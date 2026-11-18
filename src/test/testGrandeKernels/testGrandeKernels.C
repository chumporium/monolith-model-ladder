/*---------------------------------------------------------------------------*\
  testGrandeKernels.C

  Uji satuan TANPA mesh untuk dua tambahan validasi Grande dkk. 2006:
    1. isoterm NittaMultisite
    2. kernel difusi geometri BOLA di intraGrainDiffusion

  Nilai acuan dihitung terpisah (Python, bisection 300 iterasi / sistem
  linier numpy), lihat doc/validasi_fisisorpsi/PERSIAPAN_GRANDE2006.md.
  Keluar dgn kode 1 bila ada yang GAGAL.
\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "IOstreams.H"
#include "isothermModel.H"
#include "intraGrainDiffusion.H"
#include "Random.H"

using namespace Foam;

static label nFail = 0;
static label nPass = 0;

void check
(
    const std::string& name,
    const scalar got,
    const scalar ref,
    const scalar relTol
)
{
    const scalar err = mag(got - ref)/max(mag(ref), VSMALL);
    const bool ok = (err <= relTol);
    if (ok) { ++nPass; } else { ++nFail; }
    Info<< (ok ? "  LULUS  " : "  GAGAL  ") << name.c_str() << ": "
        << got << "  (acuan " << ref << ", galat rel " << err
        << ", tol " << relTol << ")" << nl;
}


int main(int argc, char *argv[])
{
    const scalar R = 8.314462618;

    // ================================================================ Nitta
    Info<< nl << "=== 1. Isoterm NittaMultisite (Grande 2006 Tabel 1, C3H6) ===" << nl;

    const scalar dH = R*Foam::log(0.0288/0.0103)/(1.0/423 - 1.0/473);

    dictionary d;
    d.add("type", word("NittaMultisite"));
    d.add("qmax", 2.857);
    d.add("a", 2.612);
    d.add("Kref", 0.0288e-3);
    d.add("Tref", 423.0);
    d.add("deltaH0", dH);

    autoPtr<isothermModel> iso = isothermModel::New("C3H6", d);

    check("q*(1500 Pa, 423 K)",  iso->qStar(1500, 423),  0.111259018697, 1e-9);
    check("q*(5630 Pa, 423 K)",  iso->qStar(5630, 423),  0.334591815531, 1e-9);
    check("q*(10000 Pa, 423 K)", iso->qStar(10000, 423), 0.498592402343, 1e-9);
    check("q*(1500 Pa, 473 K)",  iso->qStar(1500, 473),  0.0424480916594, 1e-9);
    check("q*(10000 Pa, 473 K)", iso->qStar(10000, 473), 0.235136232628, 1e-9);

    // K(473) dari deltaH0 harus kembali ke Tabel 1 (0,0103 1/kPa)
    check
    (
        "K(473 K) = dq*/dp(0)/qmax [1/kPa]",
        1e3*iso->dqStardp(0, 473)/2.857, 0.0103, 1e-9
    );

    for (const scalar p : {100.0, 1500.0, 5630.0, 20000.0, 90000.0})
    {
        for (const scalar T : {423.0, 473.0})
        {
            const scalar hp = 1e-4*p;
            const scalar fdp =
                (iso->qStar(p + hp, T) - iso->qStar(p - hp, T))/(2*hp);
            check
            (
                "dq*/dp vs beda-hingga p=" + std::to_string(int(p))
              + " T=" + std::to_string(int(T)),
                iso->dqStardp(p, T), fdp, 1e-6
            );

            const scalar hT = 1e-3;
            const scalar fdT =
                (iso->qStar(p, T + hT) - iso->qStar(p, T - hT))/(2*hT);
            check
            (
                "dq*/dT vs beda-hingga p=" + std::to_string(int(p))
              + " T=" + std::to_string(int(T)),
                iso->dqStardT(p, T), fdT, 1e-5
            );

            // Gamma = dln p/dln q
            const scalar q = iso->qStar(p, T);
            const scalar gFD = (1.0/p)/(iso->dqStardp(p, T)/q);
            check
            (
                "Gamma vs dlnp/dlnq p=" + std::to_string(int(p))
              + " T=" + std::to_string(int(T)),
                iso->thermodynamicFactor(q, T), gFD, 1e-9
            );
        }
    }

    // a = 1 harus identik Langmuir
    {
        dictionary d1(d);
        d1.set("a", 1.0);
        autoPtr<isothermModel> iso1 = isothermModel::New("C3H6", d1);
        for (const scalar p : {10.0, 1500.0, 1e5})
        {
            const scalar x = 0.0288e-3*p;
            check
            (
                "a=1 -> Langmuir p=" + std::to_string(int(p)),
                iso1->qStar(p, 423), 2.857*x/(1 + x), 1e-12
            );
        }
    }

    // ================================================================ Bola
    Info<< nl << "=== 2. Kernel difusi BOLA ===" << nl;

    // Momen pertama semi-diskret (acuan numpy) & analitik 1/15 (satuan R^2/D)
    const List<Pair<scalar>> ref
    ({
        Pair<scalar>(8,  0.06840007),
        Pair<scalar>(16, 0.06710052),
        Pair<scalar>(32, 0.06677516)
    });

    forAll(ref, k)
    {
        const label N = label(ref[k].first());
        const scalar hRad = 1.0/N;                // R = 1, D = 1
        const scalar dt = 2e-3;                    // sengaja kasar
        const scalar r = dt/(hRad*hRad);

        scalarList a(N), b(N), c(N), w(N), rhs(N), cp(N), dw(N), x(N, 0.0);
        scalar sRhs = 0;
        kineticClosures::intraGrainDiffusion::assembleSphere
        (
            N, r, a, b, c, sRhs, w
        );

        // Euler implisit dari nol, q* = 1. Momen pertama jumlah titik-kanan
        // = momen semi-diskret TEPAT (lihat catatan di intraGrainDiffusion.H).
        scalar M1 = 0, F = 0;
        label steps = 0;
        do
        {
            rhs = x;
            rhs[N - 1] += sRhs;
            F = kineticClosures::intraGrainDiffusion::solveTridiagMean
            (
                a, b, c, rhs, w, false, cp, dw, x
            );
            M1 += dt*(1 - F);
            ++steps;
        } while (1 - F > 1e-13 && steps < 10000000);

        check
        (
            "M1 bola N=" + std::to_string(N) + " vs semi-diskret",
            M1, ref[k].second(), 2e-6
        );
        check
        (
            "M1 bola N=" + std::to_string(N) + " vs analitik 1/15",
            M1, 1.0/15.0, (N == 8 ? 0.03 : N == 16 ? 0.007 : 0.002)
        );
        check("kejenuhan akhir F", F, 1.0, 1e-12);
    }

    // Konservasi satu langkah dari profil acak:
    //   Sum w_i (x_i - old_i) = sRhs*w_{N-1} (q* - x_{N-1})
    {
        const label N = 16;
        Random rnd(1234);
        scalarList a(N), b(N), c(N), w(N), rhs(N), old(N), cp(N), dw(N), x(N);
        scalar sRhs = 0;
        kineticClosures::intraGrainDiffusion::assembleSphere
        (
            N, 3.7, a, b, c, sRhs, w
        );
        forAll(old, i) { old[i] = rnd.sample01<scalar>(); }
        const scalar qs = 0.83;
        rhs = old;
        rhs[N - 1] += sRhs*qs;
        kineticClosures::intraGrainDiffusion::solveTridiagMean
        (
            a, b, c, rhs, w, false, cp, dw, x
        );
        scalar dMass = 0;
        forAll(x, i) { dMass += w[i]*(x[i] - old[i]); }
        const scalar flux = sRhs*w[N - 1]*(qs - x[N - 1]);
        check("konservasi 1 langkah (massa vs fluks permukaan)", dMass, flux, 1e-12);
    }

    Info<< nl << "RINGKASAN: " << nPass << " lulus, " << nFail << " gagal" << nl
        << endl;

    return (nFail == 0) ? 0 : 1;
}


// ************************************************************************* //
