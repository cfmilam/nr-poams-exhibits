/* Atomic ledger: inherited isotope labels and a calibrated circular pair.
   The two computations are deliberately separate: a mass-number label is
   neither a measured endpoint mass nor a physical core-radius prescription. */
(function (root) {
  'use strict';
  const calibrated = typeof module === 'object' && module.exports
    ? require('./coulomb_calibrated.js') : root.POAMSCoulomb;
  function counts(Z, A) {
    if (!Number.isInteger(Z) || !Number.isInteger(A) || Z < 1 || A < Z) {
      throw new RangeError('The count ledger requires integer A ≥ Z ≥ 1');
    }
    return {Z, A, depth:A-Z, pattern:Z, ratio:(A-Z)/Z, depthFraction:(A-Z)/A};
  }
  function pair(model, n) {
    if (!['pair','equal'].includes(model)) throw new RangeError('Choose finite or equal endpoints');
    const q=calibrated.calculate(model,n,1,2), k=calibrated.constants;
    const m1=k.lightMass, m2=model==='pair'?k.heavyMass:k.lightMass;
    const M=m1+m2, hbar=k.h/(2*Math.PI), L=n*hbar;
    const share1=m2/M, share2=m1/M;
    const r1=share1*q.radius, r2=share2*q.radius;
    // Direct mechanical sums, not definitions from the displayed shares.
    const T1=m1*(q.omega*r1)**2/2, T2=m2*(q.omega*r2)**2/2;
    const L1=m1*r1*r1*q.omega, L2=m2*r2*r2*q.omega;
    return {...q,m1,m2,M,L,share1,share2,r1,r2,T1,T2,L1,L2,
      U:-q.K/q.radius, f:q.omega/(2*Math.PI), bindingRate:q.binding/k.h,
      S1:2*Math.PI*L1,S2:2*Math.PI*L2};
  }
  function mount(elements) {
    const $=id=>document.getElementById(id);
    if (!$('atomic-ledger')) return;
    const fmt=(x,d=7)=>x.toPrecision(d).replace(/(?:\.0+|(\.\d*?)0+)(e|$)/,'$1$2');
    const option=(value,label)=>new Option(label,String(value));
    elements.forEach(e=>$('al-element').add(option(e.n,`${e.n} · ${e.name} (${e.symbol})`)));
    const body=$('al-all-counts');
    elements.forEach(e=>{
      const q=counts(e.n,e.massNumber), row=body.insertRow();
      [`${e.symbol}-${q.A}`,q.Z,q.depth,q.pattern,fmt(q.ratio,5)].forEach(value=>{
        const cell=row.insertCell();cell.textContent=String(value);
      });
    });
    function isotopeOptions() {
      const e=elements.find(e=>e.n===Number($('al-element').value));
      const values=[...new Set([e.massNumber,...e.isotopes.map(s=>Number(s.match(/-(\d+)/)[1]))])];
      $('al-isotope').replaceChildren(...values.map(A=>option(A,`${e.symbol}-${A}`)));
      renderCounts();
    }
    function renderCounts() {
      const q=counts(Number($('al-element').value),Number($('al-isotope').value));
      $('al-count-summary').textContent=`${q.depth} depth : ${q.pattern} pattern — ratio ${fmt(q.ratio)} : 1`;
      $('al-count-fractions').textContent=`${fmt(100*q.depthFraction,5)}% depth and ${fmt(100*(1-q.depthFraction),5)}% pattern, of ${q.A} assigned entries.`;
      $('al-depth-bar').style.width=`${100*q.depthFraction}%`;
      $('al-pattern-bar').style.width=`${100*(1-q.depthFraction)}%`;
      $('al-count-note').textContent=q.depth===0
        ? 'Zero surplus depth entries does not mean zero intrinsic turning or no compact structure.'
        : 'These are assigned counts, not measured regional masses, energy fractions or resultant angular momenta.';
    }
    function renderPair() {
      const q=pair($('al-mass').value,Number($('al-n').value)), k=calibrated.constants;
      $('al-n-label').textContent=String(q.n);
      $('al-radius').textContent=fmt(q.radius*1e12)+' pm';
      $('al-binding-rate').textContent=fmt(q.bindingRate/1e15)+' PHz';
      $('al-circulation').textContent=fmt(q.f/1e15)+' PHz';
      $('al-binding-units').textContent=fmt(q.binding)+' J = '+fmt(q.bindingEV)+' eV';
      $('al-virial').textContent=`T/h = +${fmt((q.T1+q.T2)/k.h/1e15)} PHz; U/h = −${fmt(-q.U/k.h/1e15)} PHz; E/h = −${fmt(q.bindingRate/1e15)} PHz.`;
      $('al-rate-ratio').textContent=`νB / fcirc = n/2 = ${q.n/2}. A binding-rate equivalent is not the circulation frequency.`;
      $('al-endpoint-label').textContent=q.model==='pair'?'Finite light/heavy calibration':'Equal-light-end diagnostic';
      [['1',q.r1,q.T1,q.L1,q.S1,q.share1],['2',q.r2,q.T2,q.L2,q.S2,q.share2]].forEach(([id,r,T,L,S,share])=>{
        $(`al-r${id}`).textContent=fmt(r*1e12)+' pm';
        $(`al-t${id}`).textContent=fmt(T/k.h/1e15)+' PHz';
        $(`al-j${id}`).textContent=fmt(L/(k.h/(2*Math.PI)))+' ℏ';
        $(`al-s${id}`).textContent=fmt(S/k.h)+' h';
        $(`al-share${id}`).textContent=fmt(share*100)+'%';
      });
      $('al-pair-total').textContent=`One pair: L = ${q.n} ℏ; angular action per full circuit = ${q.n} h; r/[ℏ/(μc)] = ${fmt(q.radiusRatio,10)}.`;
    }
    $('al-element').value='6';
    $('al-element').addEventListener('change',isotopeOptions);
    $('al-isotope').addEventListener('change',renderCounts);
    ['al-mass','al-n'].forEach(id=>$(id).addEventListener('input',renderPair));
    isotopeOptions();renderPair();
  }
  const api=Object.freeze({counts,pair,mount});
  root.POAMSAtomicLedger=api;
  if(typeof module==='object' && module.exports) module.exports=api;
})(typeof globalThis!=='undefined'?globalThis:this);
