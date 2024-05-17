1. The researcher needs to store the protein sequence to be predicted into a “.txt” file, the first line stores the name of the protein sequence, the second line is the long sequence of the protein to be predicted.
   
   ![image](https://github.com/tomatodyk/Glycation-predictor/assets/107628699/d385a2f0-2e4b-40d4-863f-a348d8f83f3d)
2. The researchers then opened the downloaded predictor software, clicked the “Upload” button, and selected the corresponding protein sequence file.
   
   ![image](https://github.com/tomatodyk/Glycation-predictor/assets/107628699/a4509a20-d256-4679-baa5-daee5d0308ee)
   ![image](https://github.com/tomatodyk/Glycation-predictor/assets/107628699/ad0ed15b-3a1c-4442-a62d-d6649d3b7b02)
3. After the file is uploaded, the program will automatically parse the protein sequence, and after parsing, it will cut the whole sequence with a window size of “31”, and if it is less than the window size, it will be supplemented with “O”. Eventually, the whole sequence will be cut into multiple sub-sequences with a window size of 31, and the trained model will be called to make predictions. The predicted content will be presented on the interface, including the number of the position where glycation occurs (which position in the whole sequence has a higher probability of lysine glycation), the predicted probability (with a threshold of 0.5, greater than 0.5 is predicted to be glycation, less than 0.5 is predicted to be non-glycation), and the predicted result (glycation or non-glycation).

   ![image](https://github.com/tomatodyk/Glycation-predictor/assets/107628699/29e6ddbc-245e-441b-9927-56c781b895ef)
4. After the prediction is done, click on the “Export” button, select the location where you want to store it, name the file, and then the file will save the final prediction in the form of excel.

   ![image](https://github.com/tomatodyk/Glycation-predictor/assets/107628699/14a9f231-8f67-4d5a-a02a-1baf5a912a64)
   ![image](https://github.com/tomatodyk/Glycation-predictor/assets/107628699/e7f86326-3b94-4c7d-b408-c103cfb1f80b)
5. The final forecast content will be saved to the excel table, find the pre-saved location, you can view the content.

   ![image](https://github.com/tomatodyk/Glycation-predictor/assets/107628699/013f8682-bb5e-4390-b685-135f054bb0f7)
   ![image](https://github.com/tomatodyk/Glycation-predictor/assets/107628699/841da00b-940e-4f34-980b-2816703f08e1)

The file “example.txt” in the instance folder can be used as a reference.








