## BRIEF

The RunSVM.m MATLAB script uses LIBSVM (Library for Support Vector Machines) to train and test whether information is available in neural population responses. The data we acquired from this novel, in vivo experiment showed a low, but statistically significant, amount of information about stimulus order (temporal encoding). Further investigation is required to validate this result. Results of this experiment are unpublished, however some interesting findings were presented by Joel Yancey at the 2015 UCLA Neuroscience Undergraduate Poster Fair. His presentation is titled "Decoding Stimulus Features from Cortical Population Response."

## POSTER ABSTRACT

The cortex encodes sensory information via the spatiotemporal activation of a large population of neurons. It is important to distinguish between “spatial information,” which refers to information encoded by which sensory afferents are activated, and “temporal information,” which refers to information encoded by the temporal pattern of the activated afferents. Speech recognition is an example of a task that requires both spatial and temporal processing in the range of tens to hundreds of milliseconds.

In this experiment, we decoded stimulus features from a population of neural responses, to characterize the spatial and temporal information available (Fig. 1). Specifically, we presented different auditory tones to awake, behaving animals, while recording extracellularly from a population of neurons in the auditory cortex. One computational model, the State-Dependent Network (SDN) model, predicts that the population of cortical neurons can encode past and present sensory events. Our goal was in part to test this hypothesis by seeing if the population response encodes not only spatial information, but information about the order of the tones as well. As such, we presented paired tones separated by 100 ms and attempted to decode the pitch of the first tone from the population response to the second tone. We used a support vector machine (SVM) algorithm to determine the amount of information the neural population had about which stimuli were presented. We found that the population encoded the features of the current stimulus well. There was a low, but statistically significant, amount of information about the previous stimulus.

## ORIGINAL POSTER

"Decoding Stimulus Features From Cortical Population Responses" (Yancey, J., Halladay, L., DeGuzman, R., Blair, T., & Buonomano, D.) Poster presented at UCLA Neuroscience Undergraduate Poster Fair in California, Los Angeles.
    
## CITATIONS
   
LIBSVM - Library For Support Vector Machines
Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for support vector machines. ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011. Software available at http://www.csie.ntu.edu.tw/~cjlin/libsvm 
