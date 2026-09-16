/* Display saved calculations only; never scales helium into another atom. */
(() => {
  'use strict';
  const campaign=window.POAMS_ALL_ELEMENTS;
  if(!campaign)return;
  const $=id=>document.getElementById(id),select=$('wa-element');
  if(!select)return;
  const set=(id,value)=>{$(id).textContent=value;};
  const fmt=(v,n=6)=>Number.isFinite(v)?Number(v).toPrecision(n):'not certified';
  const u=campaign.units,rate=u.rate_Hz/1e15,aPM=u.length_m*1e12;
  const counts=campaign.counts;
  set('ae-campaign-status',`${counts.states}/118 stationary RI-JK references · ${counts.two_basis} two-basis comparisons · ${counts.relaxed_responses} checked static responses · ${counts.relaxed_removals ?? 0} separately relaxed removal intervals. ${counts.basis_sensitive_response} responses are basis-sensitive; ${counts.direct_gradient_flagged} selected references need direct-integral refinement. Missing or unresolved quantities are not filled by interpolation.`);
  const current=()=>campaign.elements[select.value];
  const rows=(id,values)=>$(id).replaceChildren(...values.map(cells=>{
    const row=document.createElement('tr');
    for(const cell of cells){const td=document.createElement('td');td.textContent=cell;row.append(td);}
    return row;
  }));
  function frame(id,xmax,ymin,ymax,xlabel,ylabel,xmin=0){
    const canvas=$(id),w=canvas.clientWidth,h=canvas.clientHeight;
    if(!w||!h)return null;
    const dpr=Math.min(window.devicePixelRatio||1,3);
    canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);
    const c=canvas.getContext('2d');c.setTransform(dpr,0,0,dpr,0,0);
    const light=document.body.getAttribute('data-theme')==='light';
    const colors={bg:light?'#f8f6ef':'#111820',text:light?'#332d21':'#d4dce4',grid:light?'#d9d4c7':'#30404a',gold:light?'#87630b':'#e3b958',blue:light?'#087792':'#59c7e0'};
    c.fillStyle=colors.bg;c.fillRect(0,0,w,h);c.font='11px system-ui';
    const left=59,right=17,top=39,bottom=45;
    const x=v=>left+(v-xmin)/(xmax-xmin)*(w-left-right),y=v=>h-bottom-(v-ymin)/(ymax-ymin)*(h-top-bottom);
    for(let j=0;j<=4;j++){
      const v=ymin+(ymax-ymin)*j/4;c.strokeStyle=colors.grid;c.beginPath();c.moveTo(left,y(v));c.lineTo(w-right,y(v));c.stroke();
      c.fillStyle=colors.text;c.textAlign='right';c.fillText(Number(v.toPrecision(2)).toString(),left-7,y(v)+4);
      const xx=xmin+(xmax-xmin)*j/4;c.textAlign='center';c.fillText(Number(xx.toPrecision(3)).toString(),x(xx),h-bottom+18);
    }
    c.textAlign='left';c.fillText(ylabel,left,19);c.textAlign='center';c.fillText(xlabel,(left+w-right)/2,h-8);
    return {c,x,y,colors,w,h,left,right,top,bottom};
  }
  function curve(f,points,color){
    const {c}=f;c.save();c.beginPath();c.rect(f.left,f.top,f.w-f.left-f.right,f.h-f.top-f.bottom);c.clip();
    c.strokeStyle=color;c.lineWidth=2;c.beginPath();points.forEach(([x,y],i)=>i?c.lineTo(f.x(x),f.y(y)):c.moveTo(f.x(x),f.y(y)));c.stroke();c.restore();
  }
  function drawRadial(){
    const r=current();if(!r)return;
    const limit=Math.max(r.r90_pm*1.7,10),points=r.radial_profile.filter(p=>p[0]<=limit);
    const peak=Math.max(...points.map(p=>p[1]));
    const f=frame('ae-radial',limit,0,peak*1.08,'Distance from compact point (pm)','Radial probability per pm');
    if(!f)return;curve(f,points,f.colors.gold);
    curve(f,[[r.mean_pm,0],[r.mean_pm,peak*.9]],f.colors.blue);
  }
  function drawResponse(){
    const r=current();if(!r||r.response.relaxed===null)return;
    const direction=$('ae-direction').value;
    const chi=direction==='average'?r.response.relaxed:r.response.tensor[Number(direction)][Number(direction)];
    const F=Number($('ae-field').value)*1e-4;
    set('ae-field-label',(F*1e4).toFixed(2));
    const shift=-.5*chi*F*F*u.rate_Hz/1e9;
    set('ae-shift',`Δ(E/h) = ${Math.abs(shift)<1e-12?'0':fmt(shift)} GHz · χ = ${fmt(chi)}`);
    set('ae-displacement',`Mean relative displacement per mobile entry: ${fmt(-chi*F*aPM/r.z)} pm. Internal response, not apparatus translation.`);
    const max=.5*chi*1e-8*u.rate_Hz/1e9;
    const f=frame('ae-curve',1,-max*1.1,max*.08,'Perturbation F (×10⁻⁴)','Energy-rate shift (GHz)',-1);
    if(!f)return;
    const points=Array.from({length:121},(_,i)=>{const x=-1+i/60;return[x,-max*x*x];});
    curve(f,points,f.colors.gold);f.c.fillStyle=f.colors.blue;f.c.beginPath();f.c.arc(f.x(F*1e4),f.y(shift),4,0,2*Math.PI);f.c.fill();
  }
  function render(){
    const r=current();$('ae-results').hidden=!r;
    if(!r)return;
    set('ae-title',`${r.symbol} — calculated all-entry reference`);
    set('ae-method',`${r.z} mobile entries · scalar relativistic RI-JK mean field · ${r.level===3?'triple':'double'}-zeta + diffuse response space`);
    set('ae-energy',fmt(r.energy_rate_PHz)+' PHz');set('ae-radius',fmt(r.mean_pm)+' pm');set('ae-rms','RMS: '+fmt(r.rms_pm)+' pm');
    set('ae-response',fmt(r.response.relaxed));
    const change=r.comparison.response_percent;
    const sensitive=Number.isFinite(change)&&Math.abs(change)>5;
    set('ae-response-grade',r.response.relaxed===null?'Response withheld: checks did not certify it':sensitive?'Basis-sensitive reference — refinement needed':'Self-consistently relaxed, orientation average');
    const compared=Object.keys(r.comparison).length>0;
    const integralSensitive=r.direct_gradient_norm>=3e-5;
    const configurationChanged=r.calculated_angular_counts.some((n,i)=>n!==r.input_angular_counts[i]);
    set('ae-quality',`${compared?'Two basis sizes checked.':'One validated basis size available.'} ${r.response.relaxed===null?'The state is calculated; no validated relaxed-response value is supplied.':sensitive?'The response changes by more than 5% between basis sizes; do not treat it as a converged molecular input.':'A solved mean-field reference, not a correlated precision result.'} ${r.response.stability_diagnosis?r.response.stability_diagnosis+' ':''}${configurationChanged?'The dominant angular labels changed during relaxation; near-equal mixed channels can change the count without a discrete occupation transfer. Compare the supplied and resulting counts below. ':''}${integralSensitive?'The direct-integral gradient also requires refinement; stationarity here applies to the RI-JK approximation. ':''}The physical core/envelope partition remains outside this effective point-component model.`);
    set('ae-radius-note',`The distribution is normalized per mobile entry and includes every occupied radial mode. The blue marker is the mean distance; the 90% radius is ${fmt(r.r90_pm)} pm. Neither is a hard outer boundary or a compact-core radius.`);
    rows('ae-energy-table',[
      ['One-entry turning and compact coupling',fmt(r.one_body_energy*rate)],
      ['Joint direct contribution',fmt(r.direct_energy*rate)],
      ['Exchange contribution from the state rule',fmt(r.exchange_energy*rate)],
      ['Total external account',fmt(r.energy_rate_PHz)]
    ]);
    set('ae-angular',`Specified 2Mₛ = ${r.input_spin_2Ms}. Calculated ⟨S²⟩/ℏ² = ${fmt(r.spin_squared)}; scalar-representation ⟨L²⟩/ℏ² = ${fmt(r.orbital_L_squared)}. Input angular counts (s,p,d,f): ${r.input_angular_counts.join(', ')}. Dominant occupied counts after relaxation: ${r.calculated_angular_counts.join(', ')}.`);
    set('ae-koopmans',`Frozen-orbital removal estimate: ${fmt(r.koopmans_removal*rate)} PHz. This is −ε of the highest occupied orbital, not a separately relaxed first-removal interval.`);
    set('ae-removal',r.removal?.certified?`Separately relaxed removal interval: ${fmt(r.removal.interval*rate)} PHz (${fmt(r.removal.interval*u.energy_eV)} eV as a conversion). This compares the neutral reference with the lowest converged +1 state among the tested spin-removal seeds; correlation and finer physical corrections remain omitted.`:'A separately relaxed removal interval has not been certified for this reference.');
    const groups=new Map();for(const o of r.occupied_orbitals){const key=o.spin+':'+o.l;if(!groups.has(key))groups.set(key,[]);groups.get(key).push(o.energy);}
    rows('ae-occupied',[...groups].map(([key,es])=>{const [s,l]=key.split(':');return[s==='0'?'A':'B',l,es.length,fmt(Math.min(...es)*rate)+' to '+fmt(Math.max(...es)*rate)];}));
    rows('ae-basis-table',Object.entries(r.levels).map(([level,v])=>[level==='2'?'Double-zeta + diffuse':'Triple-zeta + diffuse',v.validated?fmt(v.energy*rate):'not certified',v.validated?fmt(v.mean_pm):'not certified',fmt(v.response)]));
    set('ae-basis-note',compared?`Larger-minus-smaller basis: energy ${fmt(r.comparison.energy_difference*rate)} PHz; mean radius ${fmt(r.comparison.mean_radius_percent,4)}%; response ${Number.isFinite(change)?fmt(change,4)+'%':'comparison unavailable'}. This difference is a sensitivity check, not a physical error bar.`:'A second certified basis result is not available for this state.');
    set('ae-validation',`Saved representation: ${r.nao} AO functions, ${r.nmo_per_spin} retained orbitals per spin.${r.nmo_per_spin<r.nao?' Near-dependent overlap modes below the 10⁻⁶ cutoff were removed; state coefficients are rectangular.':''} RI-JK state gradient norm: ${r.gradient_norm.toExponential(2)}. Direct four-centre gradient norm: ${r.direct_gradient_norm.toExponential(2)}. Radial count error: ${r.radial_normalization_error.toExponential(2)} entries. Direct-integral energy correction at this state: ${r.integral_energy_difference.toExponential(2)} E*. ${integralSensitive?'The direct gradient exceeds the same 3×10⁻⁵ stationarity threshold; a small total-energy correction does not by itself certify the direct-integral state. ':''}These diagnostics are retained in the downloadable records.`);
    $('ae-state-download').href='./scripts/'+r.state_archive;
    set('ae-state-download',`${r.symbol} state vectors (archive)`);
    $('ae-response-panel').hidden=r.response.relaxed===null;
    $('ae-hhe-note').hidden=r.z>2;
    set('ae-frozen-note',Number.isFinite(r.response.frozen)?`Separate frozen-potential reference: χ = ${fmt(r.response.frozen)}, C₆(${r.symbol}–${r.symbol}) = ${fmt(r.response.frozen_C6)} E*a*⁶. These use an unrelaxed response spectrum; they are not the dispersion coefficients of the relaxed static tensor or the correlated H/He calculation.`:'A frozen-potential response spectrum was not certified for this state.');
    schedule();
  }
  let queued=false;function schedule(){if(queued)return;queued=true;requestAnimationFrame(()=>{queued=false;drawRadial();drawResponse();});}
  select.addEventListener('change',render);
  for(const id of ['ae-field','ae-direction'])$(id).addEventListener('input',drawResponse);
  window.addEventListener('hashchange',render);window.addEventListener('resize',schedule);
  new ResizeObserver(schedule).observe($('ae-results'));
  new MutationObserver(schedule).observe(document.body,{attributes:true,attributeFilter:['data-theme']});
  render();
})();
