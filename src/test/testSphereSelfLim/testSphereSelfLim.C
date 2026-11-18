/*---------------------------------------------------------------------------*\
  testSphereSelfLim.C

  Uji satuan TANPA mesh untuk kernel bola D(q) (assembleSphereVar), REGISTER
  "Tahap 4 (pendahuluan)", gerbang kode:
    U1  D seragam -> koefisien identik assembleSphere          (rel < 1e-14)
    U2  konservasi: sum w (x - old) = w_{N-1} sRhs (q* - x_{N-1}) (rel < 1e-12)
        dengan D acak antar kulit (kontras 1e-3..1)
    U3  alpha -> 0: q-bar jalur D(q) = jalur D konstan          (rel < 1e-12)
  Keluar dgn kode 1 bila ada yang GAGAL.
\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "intraGrainDiffusion.H"
#include "Random.H"

using namespace Foam;
typedef kineticClosures::intraGrainDiffusion IGD;

static label nFail = 0, nPass = 0;

void check(const std::string& name, const scalar got, const scalar ref, const scalar tol)
{
    const scalar err = mag(got - ref)/max(mag(ref), VSMALL);
    const bool ok = (err <= tol);
    if (ok) { ++nPass; } else { ++nFail; }
    Info<< (ok ? "  LULUS  " : "  GAGAL  ") << name.c_str() << ": " << got
        << "  (acuan " << ref << ", galat rel " << err << ", tol " << tol << ")" << nl;
}

int main(int argc, char *argv[])
{
    const label N = 16;
    const scalar R = 4.2e-9, D0 = 1.9e-20, dt = 0.0107;
    const scalar h = R/N, dth2 = dt/(h*h);

    scalarList a(N), b(N), c(N), w(N), a0(N), b0(N), c0(N), w0(N);
    scalarList rhs(N), cp(N), dw(N), x(N), Dsh(N), old(N);
    scalar s = 0, s0 = 0;

    // ---- U1: D seragam ----
    Info<< "U1: D seragam vs assembleSphere" << nl;
    for (label i = 0; i < N; ++i) { Dsh[i] = D0; }
    IGD::assembleSphereVar(N, dth2, Dsh, D0, a, b, c, s, w);
    IGD::assembleSphere(N, D0*dth2, a0, b0, c0, s0, w0);
    scalar eMax = 0;
    for (label i = 0; i < N; ++i)
    {
        eMax = max(eMax, mag(a[i] - a0[i])/max(mag(b0[i]), VSMALL));
        eMax = max(eMax, mag(b[i] - b0[i])/mag(b0[i]));
        eMax = max(eMax, mag(c[i] - c0[i])/max(mag(b0[i]), VSMALL));
        eMax = max(eMax, mag(w[i] - w0[i])/mag(w0[i]));
    }
    eMax = max(eMax, mag(s - s0)/mag(s0));
    check("maks galat koefisien", 1.0 + eMax, 1.0, 1e-14);

    // ---- U2: konservasi, D acak (kontras 1e-3..1), 5 realisasi ----
    Info<< "U2: konservasi D acak" << nl;
    Random rnd(12345);
    for (label k = 0; k < 5; ++k)
    {
        for (label i = 0; i < N; ++i)
        {
            Dsh[i] = D0*Foam::pow(10.0, -3.0*rnd.sample01<scalar>());
            old[i] = 0.02*rnd.sample01<scalar>();
        }
        const scalar Ds = D0*Foam::pow(10.0, -3.0*rnd.sample01<scalar>());
        const scalar qs = 0.024;
        // dt besar supaya perubahan tidak kecil (uji tak trivial)
        IGD::assembleSphereVar(N, 1e4*dth2, Dsh, Ds, a, b, c, s, w);
        for (label i = 0; i < N; ++i) { rhs[i] = old[i]; }
        rhs[N - 1] += s*qs;
        IGD::solveTridiagMean(a, b, c, rhs, w, false, cp, dw, x);
        scalar dM = 0;
        for (label i = 0; i < N; ++i) { dM += w[i]*(x[i] - old[i]); }
        const scalar flux = w[N - 1]*s*(qs - x[N - 1]);
        check("realisasi " + std::to_string(k) + ": sum w dq = fluks permukaan", dM, flux, 1e-12);
    }

    // ---- U3: alpha -> 0 (faktor exp(-1e-12 theta)) == D konstan ----
    Info<< "U3: alpha->0 vs D konstan" << nl;
    for (label i = 0; i < N; ++i)
    {
        old[i] = 0.01*scalar(i)/scalar(N);
        Dsh[i] = D0*Foam::exp(-1e-12*old[i]/0.038);
    }
    const scalar qs = 0.024;
    IGD::assembleSphereVar(N, dth2, Dsh, D0*Foam::exp(-1e-12*qs/0.038), a, b, c, s, w);
    for (label i = 0; i < N; ++i) { rhs[i] = old[i]; }
    rhs[N - 1] += s*qs;
    const scalar qVar = IGD::solveTridiagMean(a, b, c, rhs, w, true, cp, dw, x);
    IGD::assembleSphere(N, D0*dth2, a0, b0, c0, s0, w0);
    for (label i = 0; i < N; ++i) { rhs[i] = old[i]; }
    rhs[N - 1] += s0*qs;
    const scalar qCon = IGD::solveTridiagMean(a0, b0, c0, rhs, w0, true, cp, dw, x);
    check("q-bar D(q) vs D konstan", qVar, qCon, 1e-12);

    Info<< nl << "LULUS " << nPass << ", GAGAL " << nFail << nl << endl;
    return (nFail > 0) ? 1 : 0;
}
