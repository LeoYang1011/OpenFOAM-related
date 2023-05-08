    // Extracting values from dict
    waveReals_ = 
        (waveDict.lookupOrDefault("waveReals", List<scalar> (1, -1.0) ));
    waveImags_ = 
        (waveDict.lookupOrDefault("waveImags", List<scalar> (1, -1.0) ));
    wavePeriods_ = 
        (waveDict.lookupOrDefault("wavePeriods", List<scalar> (1, -1.0) ));
    waveDirs_ = 
        (waveDict.lookupOrDefault("waveDirs", List<scalar> (1, -1.0) ));

    genAbs_ = (waveDict.lookupOrDefault<bool>("genAbs", false ));
    nPaddles_ = (waveDict.lookupOrDefault<label>("nPaddles", 1));
    absDir_ = (waveDict.lookupOrDefault<scalar>("absDir", 0));

    // Check if the value of nPaddles is correct for the number of columns
    if (nPaddles_ < 1)
    {
        FatalError
            << "Check nPaddles value."
            << exit(FatalError);
    }

    if ( nPaddles_ > 1 )
    {
        nPaddles_ = decreaseNPaddles( nPaddles_, patchD, dMin, dSpan );
        reduce(nPaddles_, minOp<label>());
    }

    if ( waterDepth_ <= 0.0 )
    {
        waterDepth_ = calcWL(alphaCell, cellGroup, zSpan)[0];
    }

    if ( min(wavePeriods_) <= 0.0 )
    {
        FatalError
            << "Check wave period values."
            << exit(FatalError);
    }

    nComp_ = waveReals_.size();

    if ( nComp_ != wavePeriods_.size() || nComp_ != waveDirs_.size() || 
            nComp_ == 1 )
    {
        FatalError
            << "Check the number of components for "
            << "waveHeights, wavePeriods and waveDirs."
            << exit(FatalError);
    }

    waveLengths_ = 0.0*waveReals_;

    forAll(waveLengths_, compInd)
    {
        waveLengths_[compInd] = 
            StokesIFun::waveLength (waterDepth_, wavePeriods_[compInd]);
    }
