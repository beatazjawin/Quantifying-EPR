# Quantifying-EPR
This repository stores the codes which implement the semidefinite programs constructed in the articles "Quantifying EPR: the resource theory of nonclassicality of common-cause assemblages" and "The resource theory of nonclassicality of channel assemblages" by Beata Zjawin, David Schmid, Matty J. Hoban, Ana Belén Sainz.

Given two assemblages, the codes check whether one can be converted into the other with local operations and shared randomness. The code is written in MATLAB and it requires CVX and QETLAB. 

The codes correspond to the following types of assemblages (EPR scenarios):  
LOSR_conversions.m - standard assemblage (standard EPR scenario) with |A|=|X|=2   
LOSR_conversions_channel_X2.m - channel assemblage (channel EPR scenario) with |A|=|X|=2  
LOSR_conversions_channel_X3.m - channel assemblage (channel EPR scenario) with |A|=2, |X|=3  
LOSR_conversions_BwI_X2.m - Bob-with-input assemblage (Bob-with-input EPR scenario) with |A|=|X|=2  
LOSR_conversions_BwI_X2.m - Bob-with-input assemblage (Bob-with-input EPR scenario) with |A|=2, |X|=3  
LOSR_conversions_MDI.m - measurement-device-independent assemblage (measurement-device-independent EPR scenario) with |A|=2, |X|=3

