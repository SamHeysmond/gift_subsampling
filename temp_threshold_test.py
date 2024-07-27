# Example list of p-values
# pvalues = np.array([0.01, 0.04, 0.03, 0.05, 0.20, 0.15, 0.10])



# # Set the desired FDR level (alpha)
# alpha = 0.05

# # Number of tests
# m = len(pvalues)

# # Calculate the c(m) constant for BY procedure
# c_m = np.sum(1.0 / np.arange(1, m + 1))

# # Sort p-values
# sorted_pvalues = np.sort(pvalues)

# # Calculate the BY threshold value
# threshold_value = 0
# for i in range(m):
#     by_threshold = alpha * (i + 1) / (m * c_m)
#     if sorted_pvalues[i] <= by_threshold:
#         threshold_value = by_threshold

# # Print the BY threshold value
# print("Benjamini & Yekutieli threshold value:", threshold_value)


#====

# pvalues = np.array([0.361,0.387,0.005,0.009,0.022,0.051,0.101,0.019])

#pvalues = np.array([0.001,0.008,0.039,0.041,0.042,0.060,0.074,0.205])

# Example_pval_data_200 = pandas.read_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Mo98_GIFT_200_ALL.csv")

# Example_pval_data_1000 = pandas.read_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Mo98_GIFT_1000_ALL.csv")

# pSNP5_list_200 = Example_pval_data_200['AVERAGE_PSNP5'].tolist()

# pSNP5_list_1000 = Example_pval_data_1000['AVERAGE_PSNP5'].tolist()

# pvalues = np.array([0.001,0.008,0.039,0.041,0.042,0.060,0.074,0.3424324234])


# # expecting a smaller threshold once log transformed
# pvalues_longer = np.array([0.001,0.008,0.039,0.041,0.042,0.060,0.074,0.3424324234,0.002,0.0123,0.0124,0.04314])



# # test the lengths of GIFT 1000 all and GIFT 200 ALL if different that might be the issue



# def get_bhy_thres(pvals, fdr_thres=0.05):

#     m = len(pvals)
#     m_float = float(m)
#     s_pvals = sorted(pvals)
#     s = 1.0
#     for i, p in enumerate(s_pvals):
        
#         print(f"i value: {i} /// p value: {p}")

# 		# is it >2 because used enumerate here? or is it diff to R script for other reason?
#         if i > 2:

#             s = s + 1.0/(i-1)

#             print(f"s value updated: {s}")

#         try:

#             print(f"thes_pval old: {thes_pval}")
#         except:
#              print(f"thes_pval old not made yet")


#         # critical value
#         thes_pval = ((i + 1.0) / m_float) * fdr_thres / s

#         print(f"thes_pval updated: {thes_pval}")

#         if p > thes_pval:

#             print(f"p:{p} is bigger than thes_pval: {thes_pval}")

#             break


#     return {'thes_pval':thes_pval, 'thres_i':i}




# seems to work

# import numpy as np

# def benjamini_yekutieli(pvalues, alpha):
#     pvalues = np.array(pvalues)
#     m = len(pvalues)
    
#     # Sort p-values
#     sorted_pvalues = np.sort(pvalues)
    
#     # Calculate the harmonic sum H_m
#     H_m = np.sum(1.0 / np.arange(1, m+1))
    
#     # Calculate BY thresholds
#     by_thresholds = np.arange(1, m+1) * alpha / (m * H_m)
    
#     # Determine which p-values are significant
#     significant = sorted_pvalues <= by_thresholds
    
#     return sorted_pvalues, by_thresholds, significant

# # Example p-values
# pvalues = [0.01, 0.04, 0.03, 0.002, 0.05]
# alpha = 0.05

# sorted_pvalues, by_thresholds, significant = benjamini_yekutieli(pvalues, alpha)

# for i, (p, thr, sig) in enumerate(zip(sorted_pvalues, by_thresholds, significant)):
#     print(f"P-value {p:.4f}, Threshold {thr:.4f}, Significant: {sig}")

#################################
### POTENTIALLY WORKING METHOD

# import numpy as np
# import pandas

# import statsmodels.api as sm
# import statsmodels
# from patsy import dmatrices



# def benjamini_yekutieli(pvalues, alpha):
#     pvalues = np.array(pvalues)
#     m = len(pvalues)
    
#     # Sort p-values
#     sorted_pvalues = np.sort(pvalues)
    
#     # Calculate the harmonic sum H_m
#     H_m = np.sum(1.0 / np.arange(1, m+1))
    
#     # Calculate BY thresholds
#     by_thresholds = np.arange(1, m+1) * alpha / (m * H_m)
    
#     # Determine which p-values are significant
#     significant = sorted_pvalues <= by_thresholds
    
#     # Find the largest threshold corresponding to a significant p-value
#     if np.any(significant):
#         max_significant_threshold = by_thresholds[significant].max()
#     else:
#         max_significant_threshold = None
    
#     return sorted_pvalues, by_thresholds, significant, max_significant_threshold

# # Example p-values


# #pvalues = [0.01, 0.04, 0.03, 0.002, 0.05]
# # GIFT_pval_data_200 = pandas.read_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Mo98_GIFT_200_ALL.csv")

# GIFT_pval_data_1000 = pandas.read_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Mo98_GIFT_1000_ALL.csv")

# # pvalues= GIFT_pval_data_200["AVERAGE_PSNP5"].tolist()

# pvalues= GIFT_pval_data_1000["AVERAGE_PSNP5"].tolist()

# print(f"List looks like {pvalues[0:5]}")

# alpha = 0.05

# sorted_pvalues, by_thresholds, significant, max_significant_threshold = benjamini_yekutieli(pvalues, alpha)

# # for i, (p, thr, sig) in enumerate(zip(sorted_pvalues, by_thresholds, significant)):
# #     print(f"P-value {p:.4f}, Threshold {thr:.4f}, Significant: {sig}")

# print(f"Max significant threshold: {max_significant_threshold:.4f}" if max_significant_threshold else "No significant p-values")


## ===
##STATSMODELS
#statsmodels.stats.multitest.fdrcorrection(pvalues, alpha=0.05, method='fdr_by', is_sorted=False)
# =====


##
# from statsmodels.stats.multitest import multipletests

# # Example p-values
# p_values = [0.01, 0.04, 0.03, 0.002, 0.05, 0.2, 0.06]

# # Perform the BY correction
# reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(p_values, alpha=0.05, method='fdr_by')

# # Output the results
# print("Reject null hypothesis:", reject)
# print("Corrected p-values:", pvals_corrected)
# print("Alpha for Sidak method:", alphacSidak)
# print("Alpha for Bonferroni method:", alphacBonf)

# # Find the largest threshold and corresponding unadjusted p-value
# import numpy as np

# # Get indices where the null hypothesis is rejected
# rejected_indices = np.where(reject)[0]

# if len(rejected_indices) > 0:
#     # Get the corrected p-values where null hypothesis is rejected
#     rejected_pvals_corrected = pvals_corrected[rejected_indices]
    
#     # Find the maximum corrected p-value among the rejected hypotheses
#     max_corrected_pval = np.max(rejected_pvals_corrected)
    
#     # Find the index of this maximum corrected p-value
#     max_index = np.where(pvals_corrected == max_corrected_pval)[0][0]
    
#     # Get the corresponding unadjusted p-value
#     corresponding_unadjusted_pval = p_values[max_index]
    
#     print("Largest threshold giving a significant p-value (corrected):", max_corrected_pval)
#     print("Corresponding unadjusted p-value:", corresponding_unadjusted_pval)
# else:
#     print("No hypothesis was rejected.")

#######################


# import numpy as np

# def benjamini_yekutieli(pvalues, alpha):
#     pvalues = np.array(pvalues)
#     m = len(pvalues)
    
#     # Sort p-values
#     sorted_pvalues = np.sort(pvalues)
    
#     # Calculate the harmonic sum H_m
#     H_m = np.sum(1.0 / np.arange(1, m+1))
    
#     # Calculate BY thresholds
#     by_thresholds = np.arange(1, m+1) * alpha / (m * H_m)
    
#     # Determine which p-values are significant
#     significant = sorted_pvalues <= by_thresholds
    
#     # Find the largest threshold with a significant p-value
#     if any(significant):
#         max_significant_idx = np.where(significant)[0][-1]
#         max_threshold = by_thresholds[max_significant_idx]
#         associated_pvalue = sorted_pvalues[max_significant_idx]
#     else:
#         max_threshold = None
#         associated_pvalue = None

#     return sorted_pvalues, by_thresholds, significant, max_threshold, associated_pvalue

# # Example p-values
# pvalues = [0.01, 0.04, 0.03, 0.002, 0.05, 0.2, 0.06]
# alpha = 0.05

# sorted_pvalues, by_thresholds, significant, max_threshold, associated_pvalue = benjamini_yekutieli(pvalues, alpha)

# for i, (p, thr, sig) in enumerate(zip(sorted_pvalues, by_thresholds, significant)):
#     print(f"P-value {p:.4f}, Threshold {thr:.4f}, Significant: {sig}")

# if max_threshold is not None:
#     print(f"\nLargest threshold with significant p-value: {max_threshold:.4f}")
#     print(f"Associated p-value: {associated_pvalue:.4f}")
# else:
#     print("\nNo significant p-values found.")





# import numpy as np
# import pandas
# # from math import log10

# def benjamini_yekutieli(pvals, alpha=0.05):
#     """
#     Benjamini-Yekutieli procedure to control the false discovery rate (FDR)
    
#     Parameters:
#     pvals (array-like): List or array of p-values from the hypothesis tests
#     alpha (float): Desired FDR level
    
#     Returns:
#     list: Boolean list indicating which p-values are significant
#     """
#     pvals = np.asarray(pvals)
#     n = len(pvals)
#     sorted_indices = np.argsort(pvals)
#     sorted_pvals = pvals[sorted_indices]
    
#     # Compute the BY critical values
#     harmonic_number = np.sum(1.0 / np.arange(1, n+1))
#     by_critical_values = np.arange(1, n+1) * alpha / (n * harmonic_number)
    
#     # Determine the largest k for which p(k) <= by_critical_value(k)
#     below_threshold = sorted_pvals <= by_critical_values
#     if np.any(below_threshold):
#         max_k = np.max(np.where(below_threshold)[0])
#         threshold_pval = sorted_pvals[max_k]
#         print(f"Threshold pval = {threshold_pval}")
#         transformed_bhy_thres = (-log10(threshold_pval))
#         print(f"LOGGED Threshold pval = {transformed_bhy_thres}")
#     else:
#         threshold_pval = 0
    
#     # Create the list of boolean values indicating significance
#     # significant = pvals <= threshold_pval
#     # return significant.tolist()

# # Example usage
# pvals = [0.01, 0.04, 0.03, 0.05, 0.002, 0.07, 0.899999]
# #GIFT_pval_data_200 = pandas.read_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Mo98_GIFT_200_ALL.csv")

# # GIFT_pval_data_1000 = pandas.read_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Mo98_GIFT_1000_ALL.csv")

# # #pvals= GIFT_pval_data_200["AVERAGE_PSNP5"].tolist()

# # pvals= GIFT_pval_data_1000["AVERAGE_PSNP5"].tolist()

# alpha = 0.05

# benjamini_yekutieli(pvals, alpha)

##############################################################
############### KEEP

# import argparse, random, pandas, os
# from math import sqrt, erfc, log10
# import numpy as np


# def get_bhy_thres(pvals, fdr_thres=0.05):

#     m = len(pvals)
#     m_float = float(m)
#     print(f"m value is: {m_float}")

#     s_pvals = sorted(pvals)
#     print("Sorted pvals are as follows")
#     print(s_pvals)

#     s = 1.0

#     for i, p in enumerate(s_pvals):

#         i+=1 #Added this as both programs started with i as 0 but then R incremented to start at 1 while python stayed at 0
#             # this was causing the difference in calculations but which one was right?
#             # do we want i to start off as 0 in the first iteration or 1?

#         print(f"Current i: {i} /// current p : {p} /// current s : {s}")

# 		# is it >2 because used enumerate here? or is it diff to R script for other reason?
#         # change to i>1 to match R script
#         if i > 1:
#             s = s + (1.0/(i-1)) #added brackets yellow
#             print(f"S updated to : {s}")

#         thes_pval = ((i + 1.0) / m_float) * (fdr_thres / s) # added brackets after *
#         print(f"Current thes_pval: {thes_pval}")

#         if p > thes_pval:

#             print(f"p: {p} was bigger than thes_pval: {thes_pval}")
#             break

#     return {'thes_pval':thes_pval, 'thres_i':i}

# my_result = get_bhy_thres(pvals)
# print(f"My result from GIFT: {my_result}")

############### 
##############################################################


##############################################################
############### KEEP

import pandas
from math import log10

#GIFT_pval_data_200 = pandas.read_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Mo98_GIFT_200_ALL.csv")

GIFT_pval_data_1000 = pandas.read_csv("/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Mo98_GIFT_1000_ALL.csv")

# example_pvals= GIFT_pval_data_200["AVERAGE_PSNP5"].tolist()

example_pvals= GIFT_pval_data_1000["AVERAGE_PSNP5"].tolist()

# all should be insig (for BY)
#example_pvals = [0.00356,0.01042,0.01208,0.02155,0.03329,0.11542]

# all should be insig (for BY)
#example_pvals=[0.01,0.04,0.03,0.05,0.20,0.15,0.10,0.25,0.30,0.50]

# rank 1 should be sig (for BY)
#example_pvals=[0.001,0.01,0.03,0.05,0.07,0.10,0.12,0.15,0.20,0.25]

def BY_procedure(pvals):

    # set up the constants
    alpha=0.05 # this is the fdr_threshold we want to control

    #alpha=0.01 # this is the fdr_threshold we want to control
    k = float(len(pvals))
    BY_threshold = None
    BY_threshold_rank= None

    # cumulative sum of K needed for BY procedure
    k_sum=0
    for value in range(1,int(k)+1):
        k_sum += 1/value

    # calculate different version of alpha for BY
    alpha_prime = (alpha/k_sum)
    #print(f"alpha: {alpha} /// k: {k} /// k_sum: {k_sum} /// alpha_prime: {alpha_prime}")

    # sort pvalues into ascending order
    s_pvals = sorted(pvals)
    # print("Sorted pvals are as follows")
    # print(s_pvals[0:5])

    for i, p in enumerate(s_pvals):

        i+=1 # must start at "rank 1"

        # check my values
        #print(f"Current i: {i} /// current p : {p} ///")

        # calculate adjusted value
        #BH calculation
        # adjusted_p = alpha*(i/k) 

        # BY calculation
        adjusted_p = alpha_prime*(i/k)

        # check my values
        #print(f"adjusted_p: {adjusted_p}")

        if p<adjusted_p:

            #print(f"{p} < {adjusted_p}")

            BY_threshold = adjusted_p
            BY_threshold_rank = i

            pass
        else:
            #print(f"{p} >{adjusted_p}: LOOP BROKEN!")
            break

    return (BY_threshold, BY_threshold_rank)

BY_threshold, BY_threshold_rank=BY_procedure(example_pvals)

if BY_threshold==None:
       print(f"BY_threshold is null; no values were deemed significant")
else:
    print(f"BY_threshold: {BY_threshold} /// BY_threshold_rank : {BY_threshold_rank} ///")

    transformed_bhy_thres = (-log10(BY_threshold))
    print(f"LOGGED BY_threshold: {transformed_bhy_thres} ///")

    BF_threshold = -log10(0.05/int(1000))
    print(f"LOGGED BF_threshold: {BF_threshold} ////")

