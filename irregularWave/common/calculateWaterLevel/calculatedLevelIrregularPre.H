    forAll(calculatedLevel, indWL)
    {
        calculatedLevel[indWL] = waterDepth_;

        forAll(waveReals_, lnd)
        {
            double faseTot = 
                waveKxs[lnd]*xGroup[indWL] + waveKys[lnd]*yGroup[indWL] 
              - waveOmegas[lnd]*currTime;

            calculatedLevel[indWL] += 
                timeMult
               *(
                    waveReals_[lnd]*cos(faseTot) 
                  + waveImags_[lnd]*sin(faseTot)
                );
        }
   }
