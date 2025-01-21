# WAV2AY
Tool using Fourier transform, doing sample conversion to simple frequencies

## Multiple outputs are supported:
- raw output for any machine owning an AY
- dmalist for instant Amstrad Plus use
- cpclist for special CPC replay (packed+interrupting player)


## usage:
```
wav2ay.exe <wavfile> <options>
```
## options:
```
-preamp <value>  amplification | default 1.0
-tresh  <value>  minimal energy| default 0.25 is minimum (max usable approx 15)
-replay <value>  frequency play| default 10
-high   <value>  highcut filter| default 4000
-low    <value>  lowcut filter | default is replay (min: 20Hz)
-wfreq  <value>  AY frequency  | default 1000000 (1MHz)
-nbchan <value>  nb channel    | default 3
-chans  <value>  channel used  | default 'ABC'
-noiseblocker    disable noise management
-dmalist         output optimised DMA list
-cpclist         output optimised list for CPC replay
-wavout <file>   output WAV preview
-aki    <file>   output AKI file for Arkos Track 2
-verbose
```

![Street Fighter 2 Adoken conversion](https://github.com/EdouardBERGE/wav2ay/blob/main/mimic.png)


wav2ay evolutions were suggested by Maitre Joe for the use in his last game https://amstrad-ggp.itch.io/mightysf

![Mighty Steel Fighter 2](https://github.com/EdouardBERGE/wav2ay/blob/main/mighty.png)



