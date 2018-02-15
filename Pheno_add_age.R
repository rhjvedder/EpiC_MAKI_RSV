library(foreign)
Phenotype <- read.spss("D:\\Dropbox\\Rolf\\Supplement\\Phenotype\\Database_selection_MAKI_genomics_Groningen_4_Dec_2017.sav", to.data.frame = TRUE)
Phenotype.age <- read.csv("D:\\Dropbox\\Rolf\\Supplement\\Phenotype\\Age_at_MAKI_III_nasal_sampling_17_Jan_2018.csv")
Phenotype.aged <- apply(Phenotype, 1, function(x) {
  if (x[1] %in% Phenotype.age$ï..Trial.nummer) {
    y <- Phenotype.age[grep(x[1], Phenotype.age$ï..Trial.nummer, ignore.case=T),3]
    y <- y/52
  } else {
    y <- NA
  }
  y
})
Phenotype$Leeftijd.MAKI3bezoek_Nasal_Sample_in_years <- Phenotype.aged
write.csv(Phenotype, "D:\\Dropbox\\Rolf\\Supplement\\Phenotype\\Age_at_MAKI_III_nasal_sampling_17_Jan_2018.csv", sep = "\t")
