For a CRAN submission we recommend that you fix all NOTEs, WARNINGs and ERRORs.
## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
❯ On windows-x86_64-devel (r-devel)
  checking CRAN incoming feasibility ... [14s] NOTE
  Maintainer: 'Lauren Eyler Dang <lauren.eyler@berkeley.edu>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    CVTMLE (11:96, 11:350)
    RCT (4:9, 11:268, 11:284, 11:587)
    TMLE (3:31)
    al (11:870)
    et (11:867)
    
Comments: This is a new submission. The possibly misspelled words are abbreviations that are either common or were previously defined.

❯ On windows-x86_64-devel (r-devel)
  checking examples ... [13s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
            user system elapsed
  ES.cvtmle 5.45   0.24    5.71
  wash      5.39   0.14    5.53

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
    
Comments: This is a new submission. The possibly misspelled words are abbreviations that are either common or were previously defined.

❯ On fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... [4s/13s] NOTE
  Maintainer: ‘Lauren Eyler Dang <lauren.eyler@berkeley.edu>’
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    CVTMLE (11:96, 11:350)
    RCT (4:9, 11:268, 11:284, 11:587)
    TMLE (3:31)
    al (11:870)
    et (11:867)
    
Comments: This is a new submission. The possibly misspelled words are abbreviations that are either common or were previously defined.

❯ On fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found

0 errors ✔ | 0 warnings ✔ | 6 notes ✖
