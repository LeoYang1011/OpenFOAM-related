        if (zSup[cellIndex] <= heights[cellGroup[cellIndex]-1]) 
        {// Full cell
            forAll(waveReals_, ln2)
            {
                faseTot = 
                    waveKxs[ln2]*xGroup[cellGroup[cellIndex]-1] 
                    + waveKys[ln2]*yGroup[cellGroup[cellIndex]-1] 
                    - waveOmegas[ln2]*currTime;
               
                double numerator = min(waveKs[ln2]*zExtra,710.4);
                double denumerator = min(waveKs[ln2]*waterDepth_,710.4);

                patchU[cellIndex] +=
                    sqr(cos( limAngle(waveAngles[ln2] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    (waveReals_[ln2]*cos(faseTot) + waveImags_[ln2]*sin(faseTot))*
                    waveOmegas[ln2]*cosh(numerator)/
                    sinh(denumerator)*cos(waveAngles[ln2]);

                patchV[cellIndex] +=
                    sqr(cos( limAngle(waveAngles[ln2] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    (waveReals_[ln2]*cos(faseTot) + waveImags_[ln2]*sin(faseTot))*
                    waveOmegas[ln2]*cosh(numerator)/
                    sinh(denumerator)*sin(waveAngles[ln2]);

                patchW[cellIndex] +=
                    sqr(cos( limAngle(waveAngles[ln2] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    (waveReals_[ln2]*cos(faseTot) - waveImags_[ln2]*sin(faseTot))*
                    waveOmegas[ln2]*sinh(numerator)/sinh(denumerator);
            }
        }
        else if (zInf[cellIndex] < heights[cellGroup[cellIndex]-1]) 
        {// Partial cell
            auxiliar = heights[cellGroup[cellIndex]-1] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/2.0;
            if (!extrapolation_)
            {
                auxiliar = zExtra;
            }
        
            forAll(waveReals_, ln3)
            {
                faseTot = 
                    waveKxs[ln3]*xGroup[cellGroup[cellIndex]-1] 
                    + waveKys[ln3]*yGroup[cellGroup[cellIndex]-1] 
                    - waveOmegas[ln3]*currTime;

                double numerator = min(waveKs[ln3]*auxiliar,710.4);
                double denumerator = min(waveKs[ln3]*waterDepth_,710.4);
          
                patchU[cellIndex] += 
                    sqr(cos( limAngle(waveAngles[ln3] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    auxiliarTotal*
                    (waveReals_[ln3]*cos(faseTot) + waveImags_[ln3]*sin(faseTot))*
                    waveOmegas[ln3]*cosh(numerator)/
                    sinh(denumerator)*cos(waveAngles[ln3]);

                patchV[cellIndex] += 
                    sqr(cos( limAngle(waveAngles[ln3] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    auxiliarTotal*
                    (waveReals_[ln3]*cos(faseTot) + waveImags_[ln3]*sin(faseTot))*
                    waveOmegas[ln3]*cosh(numerator)/
                    sinh(denumerator)*sin(waveAngles[ln3]);

                patchW[cellIndex] += 
                    sqr(cos( limAngle(waveAngles[ln3] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    auxiliarTotal*
                    (waveReals_[ln3]*cos(faseTot) - waveImags_[ln3]*sin(faseTot))*
                    waveOmegas[ln3]*sinh(numerator)/sinh(denumerator);
            }
        }
        else // Empty cell
        {
            patchU[cellIndex] = 0;
            patchV[cellIndex] = 0;
            patchW[cellIndex] = 0;
        }
