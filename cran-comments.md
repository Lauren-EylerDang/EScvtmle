For a CRAN submission we recommend that you fix all NOTEs, WARNINGs and ERRORs.
## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
❯ On windows-x86_64-devel (r-devel)
  checking CRAN incoming feasibility ... [15s] NOTE
  Maintainer: 'Lauren Eyler Dang <lauren.eyler@berkeley.edu>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    CVTMLE (11:96, 11:350)
    RCT (4:9, 11:268, 11:284, 11:587)
    TMLE (3:31)
    al (11:870)
    et (11:867)
    
Comment: This is a new submission. The possibly misspelled words are abbreviations that are either common or previously defined.

❯ On windows-x86_64-devel (r-devel)
  checking examples ... [20s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
            user system elapsed
  ES.cvtmle 8.56   0.39    9.08
  wash      8.41   0.09    8.59

❯ On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

❯ On ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Lauren Eyler Dang <lauren.eyler@berkeley.edu>’
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    al (11:870)
    CVTMLE (11:96, 11:350)
    et (11:867)
    RCT (4:9, 11:268, 11:284, 11:587)
    TMLE (3:31)

Comment: This is a new submission. The possibly misspelled words are abbreviations that are either common or previously defined.


❯ On ubuntu-gcc-release (r-release)
  checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
             user system elapsed
  ES.cvtmle 8.695  0.043  18.443
  wash      8.551  0.012  17.945

❯ On fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... [6s/22s] NOTE
  Maintainer: ‘Lauren Eyler Dang <lauren.eyler@berkeley.edu>’
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    CVTMLE (11:96, 11:350)
    RCT (4:9, 11:268, 11:284, 11:587)
    TMLE (3:31)
    al (11:870)
    et (11:867)
    
Comment: This is a new submission. The possibly misspelled words are abbreviations that are either common or previously defined.


❯ On fedora-clang-devel (r-devel)
  checking examples ... [20s/42s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
             user system elapsed
  ES.cvtmle 8.942  0.018  18.547
  wash      8.781  0.013  19.170

❯ On fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found

0 errors ✔ | 0 warnings ✔ | 8 notes ✖
