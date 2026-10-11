(()=>{'use strict';
const fmt=t=>{t=Math.max(0,Math.floor(t||0));return `${Math.floor(t/60)}:${String(t%60).padStart(2,'0')}`};
const audios=[...document.querySelectorAll('.poams-listen audio')];
for(const root of document.querySelectorAll('.poams-listen')){
const a=root.querySelector('audio'),pill=root.querySelector('.pl-pill'),panel=root.querySelector('.pl-player'),play=root.querySelector('.pl-play'),seek=root.querySelector('.pl-seek'),timer=root.querySelector('.pl-timer'),error=root.querySelector('.pl-error');
const update=()=>{const d=Number.isFinite(a.duration)?a.duration:0;timer.textContent=`${fmt(a.currentTime)} / ${fmt(d)}`;seek.value=a.currentTime;seek.setAttribute('aria-valuetext',`${fmt(a.currentTime)} of ${fmt(d)}`);play.textContent=a.paused?'▶':'Ⅱ';play.setAttribute('aria-label',a.paused?'Play podcast':'Pause podcast')};
pill.onclick=()=>{panel.hidden=!panel.hidden;pill.setAttribute('aria-expanded',String(!panel.hidden));if(panel.hidden)a.pause();else if(a.readyState===0||a.error){error.hidden=true;a.preload='metadata';a.load()}};
play.onclick=async()=>{try{error.hidden=true;if(a.paused){audios.filter(x=>x!==a).forEach(x=>x.pause());await a.play()}else a.pause()}catch{error.textContent='Unable to play audio. Please try again.';error.hidden=false}};
seek.oninput=()=>{a.currentTime=Number(seek.value);update()};
a.addEventListener('loadedmetadata',()=>{if(Number.isFinite(a.duration)&&a.duration>0){seek.max=a.duration;seek.disabled=false;play.disabled=false}update()});
for(const e of ['timeupdate','play','pause','ended'])a.addEventListener(e,update);
a.addEventListener('error',()=>{error.textContent='Audio could not be loaded. Close AUDIO and reopen to retry.';error.hidden=false});
}})();
