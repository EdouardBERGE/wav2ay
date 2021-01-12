
tool using Fourier transform, doing sample conversion to simple frequencies

=========== conversion from WAV to AY registers ===========

usage: wav2ay.exe <wavfile> <options>
  
options:

-preamp <value>  amplification | default 1.0
  
-tresh  <value>  minimal energy| default 0.25 is minimum (max usable approx 15)
  
-replay <value>  frequency play| default 10
  
-high   <value>  highcut filter| default 4000
  
-low    <value>  lowcut filter | default is replay
  
-wfreq  <value>  AY frequency  | default 1000000 (1MHz)
  
-nbchan <value>  nb channel    | default 3
  
-dmalist         output optimised DMA list

-wavout <file>   output WAV preview
  
-verbose

![Street Fighter 2 Adoken conversion](https://github.com/EdouardBERGE/wav2ay/blob/main/mimic.png)
