test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("getGeneModelUrls returns links to bed.gz files", {
  urls = getGeneModelUrls()
  num_urls = length(urls)
  num_matched = length(grep("bed.gz$",urls))
  expect_equal(num_urls,num_matched)
})

test_that("getGeneLengthsForPanel can accept one URL in a list and gets the right data", {
  urls=c("http://igbquickload.org/quickload/M_musculus_Dec_2011/M_musculus_Dec_2011_refGene.bed.gz")
  data=getGeneLengthsForPanel(urls)
  correct_species="M_musculus"
  returned_species=as.character(unique(data$species))
  expect_equal(correct_species,returned_species)
})

test_that("getMedianGeneLengthsForPanel returns a numeric vector with medians", {
  species = factor(c(rep("species1",101),rep("species2",100)))
  log10length=c(1:101,1:100)
  correct_medians=c(median(log10length[1:101]),
                    median(log10length[102:201]))
  names(correct_medians)=c("species1","species2")
  testdata = data.frame("log10length"=log10length,
                        "species"=species)
  testresult=getMedianGeneLengthsForPanel(testdata)
  expect_equal(testresult["species2"],correct_medians["species2"])
})

test_that("getGenomeStructureForOneSpecies returns a data frame with the right species level", {
  u = "http://igbquickload.org/quickload/M_musculus_Dec_2011/M_musculus_Dec_2011_refGene.bed.gz"
  data = getGenomeStructureForOneSpecies(u)
  correct_species="M_musculus"
  returned_species=as.character(unique(data$species))
  expect_equal(correct_species,returned_species)
})

test_that("getGenomeStructureForOneSpecies returns a data frame with column names species, chr, and size", {
  u = "http://igbquickload.org/quickload/M_musculus_Dec_2011/M_musculus_Dec_2011_refGene.bed.gz"
  data = getGenomeStructureForOneSpecies(u)
  correct_answer=c("species","chr","size")
  returned_answer=names(data)
  expect_equal(correct_answer,returned_answer)
})

test_that("getGenomeStructureForOneSpecies removes unassembled chromosome assemblies", {
  u = "http://igbquickload.org/quickload/M_musculus_Dec_2011/M_musculus_Dec_2011_refGene.bed.gz"
  data = getGenomeStructureForOneSpecies(u)
  correct_answer=0
  returned_answer=length(grep("Un",data$chr))
  expect_equal(correct_answer,returned_answer)
})

test_that("getGenomeStructureForOneSpecies chromosome names that start with Chr for Arabidopsis", {
  u = "http://igbquickload.org/quickload/A_thaliana_Jun_2009/Araport11.bed.gz"
  data = getGenomeStructureForOneSpecies(u)
  correct_answer=7
  returned_answer=length(grep("Chr",data$chr))
  expect_equal(correct_answer,returned_answer)
})

test_that("getGenomeStructuresForPanel returns 151 rows and 3 columns", {
  data = getGenomeStructuresForPanel()
  correct_answer=c(151,3)
  returned_answer=dim(data)
  expect_equal(correct_answer,returned_answer)
})

test_that("getGenomeSizesForPanel returns correct size for fruit fly genome", {
  data = getGenomeSizesForPanel()
  correct_answer=c(137624933)
  names(correct_answer)="D_melanogaster"
  returned_answer=data["D_melanogaster"]
  expect_equal(correct_answer,returned_answer)
})

test_that("getGenomeSizesForPanel returns a numeric vector", {
  data = getGenomeSizesForPanel()
  correct_answer=TRUE
  returned_answer=is.numeric(data)
  expect_equal(correct_answer,returned_answer)
})






