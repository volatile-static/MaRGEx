from seq import fid, gre3d, cpmg, larmor, inversionRecovery, noise, shimmingSweep, a2re, b0map, oSpace

defaultsequences = {
    s.mapVals['seqName']: s for s in [
        larmor.Larmor(),
        oSpace.OSpace(),
        a2re.A2RE(),
        b0map.B0Map(),
        cpmg.TSE(),
        fid.FID(),
        gre3d.GRE3D(),
        inversionRecovery.InversionRecovery(),
        noise.Noise(),
        shimmingSweep.ShimmingSweep()
    ]
}
