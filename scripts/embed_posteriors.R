library(ROL)

args = commandArgs(trailingOnly=TRUE)

tax_level <- args[1]
MAP <- as.logical(args[2])

which_measure <- "Sigma"

which_distance <- "Riemannian"

cat("Calculating distances...\n")
calc_posterior_distances(tax_level=tax_level, which_measure=which_measure, which_distance=which_distance, MAP=MAP)

cat("Embedding posterior...\n")
embed_posteriors(tax_level=tax_level, which_measure=which_measure, MAP=MAP)

cat("Plotting ordination(s)...\n")
plot_ordination(tax_level=tax_level, which_measure=which_measure, annotation="host", MAP=MAP)
plot_ordination(tax_level=tax_level, which_measure=which_measure, annotation="samplenumber", MAP=MAP)
plot_ordination(tax_level=tax_level, which_measure=which_measure, annotation="group", MAP=MAP)
