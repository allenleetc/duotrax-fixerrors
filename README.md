# DuoTrax FixErrors #

DuoTrax/Fixerrors is a working prototype based on Ctrax/Fixerrors. It is specialized for two-fly tracking.

### Requirements ###

So far tested/developed only Windows7+R2016b. Use of Matlab R2014b or later is probably best.

### Install ###

Just git clone or download the repo. In Matlab, navigate to <DuoTraxFixErrors>/fixerrors and type 'fixerrors'.

### Open Issues ###
* When trx with wingtracking are supplied, "Flip Orientations" does not do anything intelligent with the wing fields .wing_anglel, .wing_angler, etc. These remain unchanged.


### Startup ###

On startup, you will be asked to specify the following files/information:

* Specify the (tracked) movie that is to be examined/fixed.
* Specify the (Duotrax-generated) trxfile, eg 'registered_trx.mat', for the movie.
* (Optional) If you saved progress previously, you will have the option of restarting previous work at this point.
* In the same directory as your registered_trx.mat file, there must be the following DuoTrax-generated files:
    1. roidata.mat
    2. trackingdata.mat
    3. bgdata.mat
* To supply user-generated suspicious frames, add a field .susp to the trx structure in the trxfile. For now, this should be a row vector the same size as trx(i).x, consisting of 0's or 1's. Frames where .susp==1 will be considered suspicious in the UI. (A .susp vector must be added to both tracks, ie trx(1).susp and trx(2).susp must both be vectors of 0's and 1's. They can be identical however.)
* At startup, you will be prompted for "Suspiciousness Parameters" in a large dialog box with multiple inputs. These were used by the original Fixerrors to calculate and find suspicious frames. At the moment, these values can all be set to NaN.

At this stage, the main UI should appear.

### Usage/Workflows ###

This prototype supports two main workflows to find and correct errors in the tracking.

**Workflow 1: Browsing Suspicious Sequences, Fixerrors-style**

This workflow uses the Sequence Navigation and Navigation Tools panes and is similar to the original fixerrors. All "bouts" of frames where the two flies are touching are compiled into a list of Suspicious Sequences. The UI opens with the first Sequence on display. You can play this sequence using the "Play Seq" or "Play Seq Slow" buttons; the purple buttons provide shortcuts to navigate to the start, middle, and end of the sequence.

If it appears that the fly identities are swapped during the sequence, this can be fixed using the "Edit Tools":

* Swap Identities: Hit "Go" and then follow the directions to swap the fly identities **from the selected frame until the end of the movie.**
* Flip Orientation: Hit "Go", select a fly by clicking on it, then follow directions to flip the selected fly's orientation for a range of frames.

After the sequence is correct, you can press the "Correct" button to navigate to the next Suspicious Sequence. At any time, you can use "Save Progress" to temporarily save your work. Or, when you are completely done, you can use "Export Trx and Quit" to save the new/corrected trx file.

**Workflow 2: Play Movie and Correct**

A second workflow simply lets you play through the movie, correcting issues as you see fit. You can navigate the movie using the scrollbar or by using the play buttons. The ">>" or "Fast Play" buttons will attempt to play quickly during uninteresting frames and slow down during more interesting frames. (At the moment the difference isn't that big. Also, the Stop buttons can sometimes be a little "sticky" for some reason.)

While browsing the movie, if there is an issue, again, you can use the Edit Tools as above.

### Notes ###
After making all corrections, you can play the entire movie start to finish. When making corrections using Suspicious Sequences, it is possible that the flies get exchanged when they are relatively far apart. 

This is a prototype and there are some rough edges and clunky behavior. We can iterate and streamline to improve the workflow.