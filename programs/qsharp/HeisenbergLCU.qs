namespace HamiltonianSimulation.HeisenbergLCU {
    open Microsoft.Quantum.Intrinsic;
    open Std.Diagnostics;
    open Std.Math;
    open HamiltonianSimulation.Common;

    operation Run(numSites : Int, couplingJ : Double, fieldH : Double, totalTime : Double) : Unit {
        use system = Qubit[numSites];
        for i in 0 .. numSites - 1 {
            Ry((if i % 2 == 0 { 1.0 } else { -1.0 }) * PI() / 4.0, system[i]);
        }
        let (coeffs, paulis) = HeisenbergPaulis(numSites, couplingJ, fieldH);
        let (weights, terms, tags) = LcuDataFromHamiltonian(coeffs, paulis, totalTime);
        ApplyLcuBlock(weights, terms, tags, system);
        ResetAll(system);
    }
}
