Skeleton code for the top quark long exercise at CMSvDAS2020@CERN
=================================================================
See the [exercise twiki](https://indico.cern.ch/event/886923/timetable/) for more information, and setup instructions.

This code can be tested with
```bash
bambooRun --maxFiles=1 -m tt5TeV.py:HelloWorld test.yml -o test_out/test_1
```
The [``test.yml``](test.yml) file is a shorter version (for checking if the code works)
of [``tt5TeV.yml``](tt5TeV.yml), which contains a reasonably complete set of samples to get started.

The [``tt5TeV.py``](tt5TeV.py) file contains a very simple analysis module
(it makes a single plot, of the dimuon invariant mass distribution), and
a few hints to get you started with writing a more realistic one.

Have fun! (and feel free to ask questions in the exercise Discord channels)
