// fluidFinalCorrectors.H -- dijalankan sesudah pEqn dan rho = thermo.rho(), pada jalur energi partisi
// (solveFluid.H) MAUPUN implisit (solveEnergyImplicit.H). Dipisah 2026-10-09 supaya kedua jalur sama.
// Korektor species final (REGISTER PB0-V2-3f/3g, 2026-10-04). YEqn di atas memakai rho hasil iterasi
// luar SEBELUMNYA, sedangkan langkah berikut memakai rho.oldTime() = rho sesudah pEqn final. Tiap
// langkah inventori species bergeser (rho_K - rho_{K-1}) Y; pada tekanan yang berubah (evakuasi)
// selisih itu bertanda tetap dan menumpuk. Obat: pada iterasi final, species diselesaikan ULANG
// dengan (rho, phi) final dan sumber gas pori dari pembacaan yang sama dengan pEqn final.
// Bawaan: aktif bila ada sumber gas pori (V2); kasus lama tidak berubah kecuali diminta.
if
(
    finalIter
 && pimple.getOrDefault<bool>("speciesFinalCorrector", SporeFluid.set(i))
)
{
    {
        #include "YEqn.H"
    }
    if (twoWayCoupling && SporeFluid.set(i))
    {
        SmassFluid[i] -= SporeFluid[i];   // YEqn.H mereset Smass ke bagian difusif saja
    }
    fadvLag = 1;
}

#include "decayScalar.H"   // skalar modus-peluruhan (Sh/Nu aliran tak tunak); mati bila tak diminta
