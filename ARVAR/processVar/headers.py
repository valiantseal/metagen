nuclCol = "base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end"

nuclColList = nuclCol.split(":")

nuclColSel = list(nuclColList[1].split(' ')) + nuclColList[5:9]

nuclColSel = nuclColList[0:2] + nuclColList[5:9]

rshead = "REF_POS_REFNT_TOTALCOUNT_A-COUNT_T-COUNT_C-COUNT_G-COUNT_A-(+)STR_T-(+)STR_C-(+)STR_G-(+)STR_A-(-)STR_T-(-)STR_C-(-)STR_G-(-)STR_A-POS_T-POS_C-POS_G-POS_A-MISMATCH_T-MISMATCH_C-MISMATCH_G-MISMATCH_Indel1_Indel1-COUNT_Indel1-(+)STR_Indel1-(-)STR_Indel1-POS_Indel1-MISMATCH_Indel2_Indel2-COUNT_Indel2-(+)STR_Indel2-(-)STR_Indel2-POS_Indel2-MISMATCH_Indel3_Indel3-COUNT_Indel3-(+)STR_Indel3-(-)STR_Indel3-POS_Indel3-MISMATCH"

headList = rshead.split("_")

def subsetList(indList):
  subList = []
  for i in indList:
    subList.append(headList[i])
  return subList

subList = subsetList(indList = indList)

roseColms = '$2,$17,$18,$19,$20,$25,$29,$31,$35,$37,$41'.replace("$", "")
roseColmInd = roseColms.split(",")
roseColmInd = [int(i) for i in roseColmInd]
roseFinalind = [x - 1 for x in roseColmInd]

subList = subsetList(indList = roseFinalind)
