
data("leduc_minimal")

####---- Getters ----####

scpModelNames(leduc_minimal)

scpModelFormula(leduc_minimal)

dim(leduc_minimal)
dim(scpModelInput(leduc_minimal))
dim(scpModelInput(leduc_minimal, filtered = FALSE))

head(scpModelFilterNPRatio(leduc_minimal))

dim(scpModelResiduals(leduc_minimal))
dim(scpModelResiduals(leduc_minimal, filtered = FALSE))
scpModelResiduals(leduc_minimal, join = FALSE)

scpModelEffects(leduc_minimal)
dim(scpModelEffects(leduc_minimal)$Set)
dim(scpModelEffects(leduc_minimal, filtered = FALSE)$Set)
scpModelEffects(leduc_minimal, join = FALSE)[[1]]

scpModelFilterThreshold(leduc_minimal)

####---- Setter ----####

scpModelFilterThreshold(leduc_minimal) <- 2
scpModelFilterThreshold(leduc_minimal)