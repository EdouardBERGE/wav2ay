
buildsna
bankset 0

setcpc 2

org #4000
run #4000
di

xor a
ld hl,#C713
call #C000 ; init de la musique

refait
ld b,#F5
.vbl in a,(c) : rra : jr nc,.vbl
.novbl in a,(c) : rra : jr c,.novbl

zika call #C003 ; jouer la musique par defaut

call scankeyboard
ld a,(KEY_ESC_BYTE) : and KEY_ESC_BIT : jr z,refait

ld hl,zesample
call sample_init ; sinon on appuie sur ESC on fait une demande de sample :)

jr refait

include 'keyboard.asm'

; format 
; 00 ; touche pas au canal
; 01 ; juste le volume
; 10 ; volume+frequence basse
; 11 ; tout
; x3 + 0/1 si registre mixer a envoyer

;******************************
          sample_init
;******************************
ld a,(iteration+1) : or a : ret nz ; si c'est pas 0 on joue déjà un sample
ld a,(hl) : inc hl : ld (iteration+1),a : ld (sample_current_offset+1),hl

; initialiser l'interruption pour le sample, ou plus simple
; comme ici on va juste claquer l'adresse de la routine
; de lecture à la place de l'adresse du player de musique
ld hl,sample_play : ld (zika+1),hl

; on enregistre le contexte de l'AY car on était en train de jouer de la musique
ld b,#F4 : exx : ld bc,#F680 : ld de,#92C0 : exx
unzero ld e,0
ld hl,music_context
enregistre_contexte
ld a,e : cp 12 : jr z,contexte_ok : call read_psg : inc e : jr enregistre_contexte
contexte_ok

; Enfin on initialise notre "sample"
; volume A/B/C a zero
ld c,8
volume_zero
ld hl,unzero+1
call routine_SendPSG
inc c : ld a,c : cp 11
jr nz,volume_zero
; mixer init tous canaux ouverts
ld c,7 : out (c),c : exx : out (c),e
defb #ED,#71 ; out (c),0
exx
ld c,%111000 : out (c),c : exx : out (c),c
defb #ED,#71 ; out (c),0
exx

xor a : ld (sample_play),a ; activer la routine
ret

music_context defs 12

read_psg out (c),a : exx : out (c),e
defb #ED,#71 ; out (c),0
inc b : out (c),d : dec b
ld a,#40 : out (c),a : exx : ini : ld bc,#F782 : out (c),c : ld b,#F4
ret

;******************************
         sample_play
;******************************
; tant qu'on ne passe pas par sample_init, la routine ne peut rien faire car elle sait qu'elle n'est pas initialisée
ret ; poké en RET dans ce POC mais l'idée c'est que la routine réactive la musique quand le sample est terminé :)

sample_current_offset ld hl,#1234

ld b,#F4 : exx : ld bc,#F680 : ld e,#C0 : exx
ld a,(hl) : inc hl ; qu est-ce qu on fait?
or a : jr z,NoMixer

add a : jr nc,channelAfreq
; volumeA
ld c,8 :call routine_SendPSG
channelAfreq
add a : jr nc,channelAlow
ld c,1 :call routine_SendPSG
channelAlow
ld c,0 :call routine_SendPSG

add a : jr nc,channelBfreq
; volumeB
ld c,9 :call routine_SendPSG
channelBfreq
add a : jr nc,channelBlow
ld c,3 :call routine_SendPSG
channelBlow
ld c,2 :call routine_SendPSG

add a : jr nc,channelCfreq
; volumeC
ld c,10 :call routine_SendPSG
channelCfreq
add a : jr nc,channelClow
ld c,5 :call routine_SendPSG
channelClow
ld c,4 :call routine_SendPSG

add a : jr nc,NoMixer
ld c,7 :call routine_SendPSG

NoMixer ld (sample_current_offset+1),hl

iteration ld a,0 : dec a : ld (iteration+1),a : ret nz

; code de fin doit nettoyer, relancer la musique, etc.
; ici, on se contente de remettre l'adresse du player de musique :)
ld hl,#C003 : ld (zika+1),hl

; pour finir, on recupere le contexte de l'AY tel qu'il était
; avant qu'on demande à jouer le sample
ld hl,music_context
xor a
restore_music_context ld c,a : call routine_SendPSG : inc a : cp 12 : jr nz,restore_music_context

; et on empeche a nouveau l'appel de sample_play sans qu'on passe par l'init d'un autre sample
ld a,#C9 : ld (sample_play),a ; RADICAL
ret

; A appeler avec A en valeur pour le registre
; inchange par la routine
routine_SendPSG
out (c),c : exx : out (c),e
defb #ED,#71 ; out (c),0
exx
inc b : outi : exx : out (c),c
defb #ED,#71 ; out (c),0
exx
ret

;print 'taille routine=',$-sample_init

zesample include 'sample.cpc'

org #C000
incbin 'ultrasyd_C713.bin'

;save 'replay.bin',#4000,$-#4000,DSK,'replay.dsk'


