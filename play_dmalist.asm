; assembly file supposed to be compiled with RASM
;
; this will generate a snapshot playing one DMA sample
;

Macro GenericDisableRom
and %11
or %10001100
ld b,#7F
out (c),a
Mend

; E : lower rom number
Macro GenericSelectLowerRom Adresse,Asic,VideoMode

if {Asic}==ON
rmr2_page=%11000
else
rmr2_page=0
endif

switch {adresse}
case #0000
; rmr2 déjà mis comme il faut
break
case #4000
if {Asic}==ON
print 'impossible de mettre la page en #4000 avec l\'asic actif'
stop
endif
rmr2_page=%01000
break
case #8000
if {Asic}==ON
print 'impossible de mettre la page en #4000 avec l\'asic actif'
stop
endif
rmr2_page=%10000
break
default
print 'impossible de selection l\'adresse ',{HEX}{adresse},' pour une ROM basse'
stop
endswitch

if {VideoMode}>7
print 'Mode video doit etre entre 0 et 3 uniquement'
stop
elseif {VideoMode}<0
print 'Mode video doit etre entre 0 et 3 uniquement'
stop
endif

ld a,rmr2_page+%10100000
add e
ld b,#7F
out (c),a
ld a,%10001000+{VideoMode} ; up deco, low co, mode 0
out (c),a
mend





ON  EQU 1
OFF EQU 0

macro Asic switch
if {switch}
ld bc,#7FB8
out (c),c
else
ld bc,#7FA0
out (c),c
endif
mend

;**********************************************************
;**********************************************************
;              SNAPSHOT configuration
;**********************************************************
;**********************************************************

buildsna
snaset cpc_type,4
snaset GA_ROMCFG,#8C
snaset GA_RAMCFG,#C0
snaset Z80_SP,#38
snaset Z80_IM,#1

bankset 0      ; initialisation is always starting in bank 0
run rom_init
org 0       ; very first code is always starting in adress 0

; sequence to unlock asic
unlockdata defb #ff,#00,#ff,#77,#b3,#51,#a8,#d4,#62,#39,#9c,#46,#2b,#15,#8a,#cd,#ee

; in case we need an interruption in the initialisation, we must set an interrupt vector in #38
org #38
ei
ret

rom_init

;; unlock ASIC so we can access ASIC registers (Kevin Thacker)
ld b,#bc
ld hl,unlockdata
ld e,17
.loop
inc b
outi
dec e
jr nz,.loop

Asic On
ld hl,#000 : ld (#6400),hl
ld hl,#555 : ld (#6402),hl
ld hl,#0F0 : ld (#6420),hl

megaloop

;******************************
;***** LAUNCH DMA LIST ********
;******************************
ld hl,sample
ld (#6C04),hl
ld hl,#6C06
ld (hl),0 ; registre de pause au minimum utilisable
ld hl,#6C0F
ld (hl),%10 ; DMA 1

jr $

ld hl,50*8
ld b,#F5
.novbl in a,(c) : rra : jr c,.novbl
.vbl in a,(c) : rra : jr nc,.vbl
ld a,e
ld (#6400),a

dec hl
ld a,h
or l
jr nz,.novbl

jr megaloop

;**********************************
; DMA sample MUST be aligned!!!
;**********************************
align 2
sample
defw #0700+%111000   ; ALL channels ON
defw 0,8,0,9,0,10    ; ALL volumes to ZERO
include 'sample.dma' ; generated DMA list    <= PUT YOUR DMA LIST FILE HERE (text file)
defw #0700+%111111   ; channel OFF
defw #4020           ; END of DMA list


