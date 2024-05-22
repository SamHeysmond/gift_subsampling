# CSV testing 
import pandas as pd


word = "csv"

with open(r"dummy_"+word+".csv", 'r') as fp:
        n_abs_theta_lines = len(fp.readlines())
        print('Total Number of lines:', n_abs_theta_lines)
print("end of script")

dataFrame = pd.read_csv("dummy_"+word+".csv")
#data1 = pd.DataFrame(columns=["CHROM","POS","PVAL"])


'''
print("\nBefore sorting = \n", dataFrame)
dataFrame.sort_values("PVAL", axis=0, ascending=True,inplace=True, na_position='first')
print("\nSorted CSV file (according to PVAL) = \n", dataFrame)
'''

print("\nLength of data frame (rows) is: \n", len(dataFrame))


new_row=pd.Series({"CHROM":1,"POS":123,"PVAL":0.001})
dataFrame =pd.concat([dataFrame, new_row.to_frame().T], ignore_index=True)

print("\nLength of data frame (rows) is: \n", len(dataFrame))
dataFrame.to_csv("dummy_"+word+".csv", index=False)

print("End of script")