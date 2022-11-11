summary <- c(23515773,38916665,37246177,28816302,40647210,49785528,43148740,47210606,46163562,44519615,46479321,45747060,47844829,47734092,54509111,46543616)

mean_depth <- mean(summary)

ratio <- mean_depth/summary


write(ratio,paste0("~/KAS-METTL/",group,"/05_bedtools/bedGraph/bam_ratio.txt"),sep=" ",ncolumns=8)
