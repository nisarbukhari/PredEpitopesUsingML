==================================================
Follow the below mentioned instructions:
1. Run FeatureExtarction.R 
2. This will generate the dataset like dataSet_SARS-CoV-2.CSV
3. In "dataSet_SARS-CoV-2.CSV" file add one column in the begining with heading as "Class" and paste the contents for this column from  the "Class" column of "Peptide with Label.CSV" file.  
4. Run HyDecisionTree.R
5. RunHyRandomForest.R
6. Run HyNeuralNetwork.R

For 4,5 and 6 above, you will get the follwoing in Current working directory
a)Model Evaluation results with evaluation metrics in Evaluation-Result.csv
b)Actual and Predicted labels in ActualPredicted-Result.csv
c)Confusion matric results in cm.csv
d)Over all CM results in cm-overall.csv
e)Results by Class as cm-class.csv
f)Cross validation results in CV.csv
g)Saved model as Model.RData
==================================================
