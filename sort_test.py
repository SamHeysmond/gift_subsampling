import pandas

alpha1=""
def test1(alpha1variable):
    alpha1="djdfdnsf"
    return alpha1

alpha1=test1(alpha1)

T20_absolute_theta_output=open("output_files/999_T20_absolute_theta.csv","w")
# write in header CHROM,POS,PVAL
T20_absolute_theta_output.write("CHROM,POS,PVAL"+"\n")
T20_absolute_theta_output.close()

# check the T20 files to see how many lines they have

# old code below
#with open(r"output_files/"+args.id+"_T20_pSNP4.csv", 'r') as fp:
    #n_abs_theta_lines = len(fp.readlines())
absolute_theta_values=[123,41,5,14,6,523423,65,245,6543,34,2135,654,2346,2344556,2344325,64654,4524546,2343435662,23432123,131334,123,22,551,54234,23]
position_values=[1233,44411,553,443,856,35611,928,200,1152,3424,999,565,34,1123,5435,2344,20,40,60,80,100,120,140,180,170]
chrom_values=[2,1,1,5,4,3,2,1,2,3,6,4,3,6,4,4,1,2,3,4,5,6,7,8,2]
# ===============================================================
# ABSOLUTE_THETA SECTION
# =============================================================
# if theres room for more SNPs (i.e. 20 limit not reached) 
# then add the current absolute theta value
# <20 since pandas doesnt count the header?
for index in range(0,len(absolute_theta_values)):
    absolute_relative_theta = absolute_theta_values[index]
    print("current absolute_relative_theta =", absolute_relative_theta)
    dataFrame_absolute_theta = pandas.read_csv("output_files/999_T20_absolute_theta.csv")
    n_abs_theta_lines= len(dataFrame_absolute_theta)
    if n_abs_theta_lines<20:
        # add in a line for the current line of data in the calculation
        new_row=pandas.Series({"CHROM":chrom_values[index],"POS":position_values[index],"PVAL":absolute_theta_values[index]})
        dataFrame_absolute_theta =pandas.concat([dataFrame_absolute_theta, new_row.to_frame().T], ignore_index=True)

        # re-sort the csv so lowest pval is at the top and biggest pval is last row
        dataFrame_absolute_theta.sort_values(by=["PVAL"], axis=0, ascending=False,inplace=True, na_position='first')
        #dataFrame_absolute_theta.reset_index()
    # if no room in list -> check if abs theta is bigger than current lowest abs theta in t20 list
    # by checking if its larger than last pval in file (should be the highest)
    else:

        # source highest pval from dataframe (should be last index and pval col)
        lowest_T20_abs_theta=dataFrame_absolute_theta.iloc[-1]["PVAL"]
        print("Lowest current T20 theta: ", lowest_T20_abs_theta)

        if float(absolute_relative_theta)> float(lowest_T20_abs_theta):
            print(absolute_relative_theta," // is bigger than // ", lowest_T20_abs_theta)
            # replace smallest value in list (last item in the list)
            # first remove the smallest abs theta (last item in the dictionary)
            dataFrame_absolute_theta=dataFrame_absolute_theta.drop(dataFrame_absolute_theta.index[-1])
            
            # add in a line for the current line of data in the calculation
            new_row=pandas.Series({"CHROM":chrom_values[index],"POS":position_values[index],"PVAL":absolute_theta_values[index]})
            dataFrame_absolute_theta =pandas.concat([dataFrame_absolute_theta, new_row.to_frame().T], ignore_index=True)

            # re-sort the csv so highest abs theta is at the top and lowest is at the bottom
            dataFrame_absolute_theta.sort_values(by=["PVAL"], axis=0, ascending=False,inplace=True, na_position='first')
            #dataFrame_absolute_theta.reset_index()
    dataFrame_absolute_theta.reset_index()
    dataFrame_absolute_theta.to_csv("output_files/999_T20_absolute_theta.csv", index=False)
    print(dataFrame_absolute_theta)
dataFrame_absolute_theta.reset_index()
print(dataFrame_absolute_theta)
print("number of rows is :", len(dataFrame_absolute_theta))

x = 15
y = 20
z = 21
if x<=z<=y:
    print("bingo")

biggest_two=dataFrame_absolute_theta.nlargest(2,"PVAL")
print("alpha1:", alpha1)
print("Biggest two DF: \n", biggest_two)