
speciesBact<-list.files('/home/ubuntu/bactopia/datasets/species-specific')

currentSpecies<-list.files('./test_gtdbtk/taxa')

presentSpecies<-currentSpecies[currentSpecies%in%speciesBact]

print(length(currentSpecies))
print(length(presentSpecies))