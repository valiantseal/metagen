# overlap
overlapSamplesDf = read.csv("result_tables/Ludy_metaAmpIvar_overlapSnv.csv")
overlapSamples = unique(overlapSamplesDf$Sample)

# metaseq
metaseq1 = read.csv("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-16/input.csv")
metaseq2 = read.csv("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-24/input.csv")
metaComb = unique(c(metaseq1$sample, metaseq2$sample))

#ampSeq 
amp1 = read.csv("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/input.csv")
amp2 = read.csv("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_ampseq_2023-06-16/input.csv")
ampComb = unique(c(amp1$sample, amp2$sample))

