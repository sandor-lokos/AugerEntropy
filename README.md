# AugerEntropy

## Workflow: SaveMult.C --> FitMultNBD.C --> CombinedEntropy.C --> PlotEntropies.C

## Files:

	- `SaveMult.C` : reads `CONEX` simulation files and save the `N0` multiplicity histogram from them into a rootfile. It also produce the auxiliary rootfile for the check of dependence on the number of entries.
	- `FitMultNBD.C` : Calculate the entropy of the multiplicity of the primaries separately and save it to a rootfile. The calculation is done directly on the histogram and based on a fit with a quadrupole of NBDs. This is a test to see how different the entropy of the sum and the sum of the entropy of the multiplicities.
	- `CombinedEntropy.C` : combine the multiplicity distributions with the weights or fractions published in [arXiv:2507.19326](https://arxiv.org/abs/2507.19326), then calculate the entropy and the average of the combined multiplicity distributions and save it into a rootfile.
	- `PlotEntropies.C` : plot the entropy as the function of ln<n> and its statistical and systematic uncertainties.
	
## Further details:

	- to be written