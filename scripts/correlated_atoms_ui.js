/* Saved correlated references. The 118-element mean-field layer is preserved. */
(() => {
  'use strict';
  const data=window.POAMS_CORRELATED_ATOMS;
  const $=id=>document.getElementById(id),select=$('wa-element');
  if(!data||!select||!$('ca-results'))return;
  const byZ={6:'C',7:'N',8:'O'},fmt=(x,n=6)=>Number.isFinite(x)?x.toPrecision(n):'not checked';
  const set=(id,v)=>{$(id).textContent=v;};
  for(const option of select.options){if(byZ[option.value])option.text=option.text.replace('calculated mean-field reference','correlated upgrade + mean-field reference');}
  function direction(){
    const atom=data.atoms[byZ[select.value]];if(!atom)return;
    const response=atom.selected.response.ccsdt,choice=$('ca-direction').value;
    const chi=choice==='average'?response.average:response[choice].value;
    const shift=-.5*chi*1e-8*data.units.rate_Hz/1e9;
    set('ca-direction-readout',`χ = ${fmt(chi)}. At a small directional perturbation F = 10⁻⁴, the local quadratic shift Δ(E/h) is ${fmt(shift)} GHz. This is an internal response, not centre-of-mass motion.`);
  }
  function render(){
    const atom=data.atoms[byZ[select.value]],show=Boolean(atom?.selected?.response?.ccsdt);
    $('ca-results').hidden=!show;if(!show)return;
    const s=atom.selected;
    $('ca-state-download').href='./scripts/'+atom.state_archive;
    set('ca-state-download',`${atom.symbol} states and amplitudes`);
    const packages=$('ca-state-packages');packages.replaceChildren(document.createTextNode('State packages for every archived basis: '));
    atom.state_archives.forEach((name,i)=>{
      if(i)packages.append(document.createTextNode(' · '));
      const link=document.createElement('a');link.href='./scripts/'+name;
      const layer=name.match(/-a(\d+)\.zip$/);link.textContent=layer?`${atom.symbol}, diffuse layer ${layer[1]}`:`${atom.symbol}, all archived bases`;packages.append(link);
    });
    set('wa-status','Correlated C/N/O accuracy upgrade + separately retained mean-field reference');
    set('ca-title',`${atom.symbol} — correlated energy, envelope and response`);
    set('ca-removal',fmt(s.removal_rate_PHz)+' PHz');
    set('ca-removal-conversion',fmt(s.removal_eV)+' eV as a conversion');
    set('ca-mean',fmt(s.mean_pm)+' pm');set('ca-rms','RMS: '+fmt(s.rms_pm)+' pm · CCSD density');
    set('ca-response',fmt(s.response.ccsdt.average));
    const convergence=atom.response_convergence;
    const grade=convergence.passed?'The average and both angular branches pass the initial 1% cardinal and diffuse response checks.':convergence.average_pass?'The average passes the initial 1% checks, but an angular branch remains basis-sensitive; the full directional response is not yet cleared.':'The response still needs basis refinement.';
    set('ca-grade',`${s.basis}; specified ${s.term} term. ${grade} Removal energies retain separate basis sensitivity. This is a correlated benchmark, not helium-level precision or a complete physical-core calculation.`);
    $('ca-basis-table').replaceChildren(...atom.levels.map(r=>{
      const row=document.createElement('tr');
      for(const value of [r.basis,fmt(r.removal_rate_PHz),fmt(r.mean_pm),fmt(r.response.ccsdt?.average)]){
        const cell=document.createElement('td');cell.textContent=value;row.append(cell);
      }return row;
    }));
    const b=atom.basis_check,d=atom.diffuse_check;
    set('ca-basis',b?`Last cardinal-basis change: average χ ${fmt(b.response_percent,4)}%; parallel ${fmt(b.response_components_percent.parallel,4)}%, perpendicular ${fmt(b.response_components_percent.perpendicular,4)}%; mean radius ${fmt(b.mean_radius_percent,4)}%; removal interval ${fmt(b.removal_eV_difference,4)} eV. Passing the response check does not establish a basis-limit removal energy.`:'Last cardinal-basis comparison unavailable.');
    set('ca-diffuse',d?`Adding diffuse layer ${s.augmentation} at the selected basis changes average χ by ${fmt(d.response_percent,4)}% (parallel ${fmt(d.response_components_percent.parallel,4)}%, perpendicular ${fmt(d.response_components_percent.perpendicular,4)}%). Response and energy have separate convergence checks.`:'Successive-diffuse-layer comparison unavailable.');
    const m=atom.benchmarks.removal,t=atom.benchmarks.response;
    set('ca-measured',`Independent measured removal interval (NIST): ${fmt(m.value_eV,9)} eV. Calculated-minus-measured: ${fmt(m.residual_eV,4)} eV. The measured ground J level and this nonrelativistic term are not identical physical models; the residual includes omitted physics and numerical approximation.`);
    set('ca-theory',`Independent nonrelativistic response calculation (Das & Thakkar, 1998): χ ≈ ${t.value}. That is a theoretical benchmark, not an experimental measurement; it was not fitted.`);
    const f=atom.small_basis_correlation_check;
    set('ca-fci',f?`Independent full-configuration check in ${f.basis}: CCSD(T)-minus-full-configuration removal interval = ${fmt(f.ccsdt_removal_minus_fci_eV,4)} eV. This tests correlation truncation in that smaller space; it is not applied as a correction here.`:'Small-space configuration check unavailable.');
    set('ca-controls',`H/He implementation controls are retained in the download. ${data.validation.records} state/field records passed normalization, stationary-state, response parity and step-size checks; ${data.validation.fci_checks} C/N/O full-configuration controls are included. Published exact H and correlated He-4 benchmarks are unchanged.`);
    direction();
  }
  select.addEventListener('change',render);$('ca-direction').addEventListener('change',direction);render();
})();
