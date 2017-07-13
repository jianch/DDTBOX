# DDTBOX
Decision Decoding Toolbox

Copyright (c) 2013--2017 Stefan Bode and contributors.

Unless otherwise specified, code is distributed under the GNU Public License (GPL) version 2, and documentation under a Creative Commons Attribution-Share-Alike 4.0 International license.

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.

We hope that you find the software and documentation useful.
If you publish an analysis using the toolbox, we ask that you cite us. 
For now, please cite us via GitHub, but we are currently preparing on a publication to submit to peer review. 
A sample citation would be:

Bode, S., Feuerriegel, D., Bennett, D., & Alday, P. (2017). The Decision Decoding ToolBOX (DDTBOX) -- A multivariate pattern analysis toolbox for event-related potentials. DOI: 10.5281/zenodo.48143

[![DOI](https://zenodo.org/badge/20756/DDTBOX/DDTBOX.svg)](https://zenodo.org/badge/latestdoi/20756/DDTBOX/DDTBOX)

If using a version from the master or development branch, please include the [commit hash](http://codetunnel.io/merge-vs-rebase-part-1-what-is-a-commit-hash/) in your citation so that others can identify the exact version of the code used for analyses.

# External Dependencies

The code in the toolbox depends on the functionality supplied by [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) and its specialised and often faster cousin, [LIBLINEAR](https://www.csie.ntu.edu.tw/~cjlin/liblinear/).
These libraries are included with DDTBOX from v1.0 onwards. In some cases you may need to configure MATLAB to use these external dependencies.
Please see their respective documentation for more information. 
