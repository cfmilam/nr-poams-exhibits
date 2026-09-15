/* Measured-alpha, two-ended circular scale calculator. No adjustable alpha.
   This computes the stated leading model, not molecular structure or thrust. */
(function (root) {
  'use strict';
  const constants = Object.freeze({alpha: 0.0072973525643, h: 6.62607015e-34,
    c: 299792458, eV: 1.602176634e-19, lightMass: 9.1093837139e-31,
    heavyMass: 1.67262192595e-27});
  const masses = Object.freeze({
    light: constants.lightMass,
    pair: constants.lightMass * constants.heavyMass / (constants.lightMass + constants.heavyMass),
    equal: constants.lightMass / 2
  });
  function calculate(model, n, lower, upper) {
    if (!Object.prototype.hasOwnProperty.call(masses, model)) throw new RangeError('Unknown mass calibration');
    if (!Number.isInteger(n) || n < 1 || n > 6) throw new RangeError('Radial mode must be 1–6');
    if (!Number.isInteger(lower) || !Number.isInteger(upper) || lower < 1 || lower >= upper || upper > 8) {
      throw new RangeError('Transition requires 1 ≤ lower < upper ≤ 8');
    }
    const {alpha, h, c, eV} = constants;
    const hbar = h / (2 * Math.PI), mu = masses[model], K = alpha * hbar * c;
    const referenceLength = hbar / (mu * c), referencePeriod = h / (mu * c*c);
    const radius1 = referenceLength / alpha, binding1 = mu * c*c * alpha*alpha / 2;
    const radius = n*n * radius1, binding = binding1 / (n*n);
    const omega = mu*c*c*alpha*alpha / (hbar*n*n*n), period = 2*Math.PI/omega;
    const transitionJ = binding1 * (1/(lower*lower) - 1/(upper*upper));
    const frequency = transitionJ / h, wavelength = c / frequency;
    return {model, n, lower, upper, mu, K, referenceLength, referencePeriod,
      radius1, binding1, radius, binding, bindingEV: binding/eV, omega, period,
      radiusRatio: radius/referenceLength, periodRatio: period/referencePeriod,
      transitionJ, frequency, wavelength, Rydberg: binding1/(h*c)};
  }
  const api = Object.freeze({constants, masses, calculate});
  root.POAMSCoulomb = api;
  if (typeof module === 'object' && module.exports) module.exports = api;
  if (typeof document === 'undefined') return;
  const $ = id => document.getElementById(id);
  if (!$('cal-mass')) return;
  const nice = (x, digits=7) => x.toPrecision(digits).replace(/(?:\.0+|(\.\d*?)0+)(e|$)/,'$1$2');
  function plot(q) {
    const bars = [], levels = [];
    for (let n=1;n<=6;n++) {
      const y=37+(n-1)*29, width=240*n*n/36;
      bars.push(`<text x="12" y="${y+5}" class="plot-label">n=${n}</text><rect x="55" y="${y-8}" width="${width}" height="14" rx="3" class="${n===q.n?'plot-active':'plot-bar'}"/><text x="${62+width}" y="${y+4}" class="plot-label">${n*n} r₁</text>`);
      const ly=20+175/(n*n);
      const labelY=195-(n-1)*29;
      levels.push(`<line x1="50" x2="280" y1="${ly}" y2="${ly}" class="${n===q.n?'plot-active-line':'plot-line'}"/><line x1="281" x2="323" y1="${ly}" y2="${labelY}" class="plot-guide"/><text x="330" y="${labelY+4}" class="plot-label">n=${n}</text>`);
    }
    $('radius-chart').innerHTML='<title>Relative mode radii, proportional to n squared</title>'+bars.join('');
    $('energy-chart').innerHTML='<title>Gross bound-state energies, minus binding divided by n squared</title><line x1="50" x2="280" y1="20" y2="20" class="plot-zero"/><text x="285" y="24" class="plot-label">0</text>'+levels.join('');
  }
  function update() {
    const lower=Number($('cal-lower').value);
    for (const option of $('cal-upper').options) option.disabled=Number(option.value)<=lower;
    if (Number($('cal-upper').value)<=lower) $('cal-upper').value=String(lower+1);
    const q=calculate($('cal-mass').value,Number($('cal-n').value),lower,Number($('cal-upper').value));
    $('cal-n-value').textContent=String(q.n);
    $('cal-radius').textContent=nice(q.radius*1e12)+' pm';
    $('cal-binding').textContent=nice(q.bindingEV)+' eV';
    $('cal-period').textContent=nice(q.period)+' s';
    $('cal-line').textContent=nice(q.wavelength*1e9)+' nm';
    $('cal-transition').textContent=`${q.upper} → ${q.lower}, vacuum; gross model`;
    $('cal-ratios').textContent=`rₙ / [ħ/(μc)] = ${nice(q.radiusRatio,9)}; Tₙ / [h/(μc²)] = ${nice(q.periodRatio,9)}. α is unchanged.`;
    $('cal-caption').textContent=q.model==='light'
      ? 'Light-mass limit: μ equals the measured light mass. This is a scale calculation, not a picture of a particle orbit.'
      : q.model==='pair'
        ? 'Finite two-ended pair: both calibrated masses participate through their reduced mass. Neither end is fixed at the barycentre.'
        : 'Equal-end diagnostic: both endpoints have the calibrated light mass. This illustrates reduced-mass scaling, not a complete species model.';
    plot(q);
  }
  for (const id of ['cal-mass','cal-n','cal-lower','cal-upper']) $(id).addEventListener('input',update);
  $('cal-reset').addEventListener('click',function () {
    $('cal-mass').value='pair'; $('cal-n').value='1'; $('cal-lower').value='2'; $('cal-upper').value='3'; update();
  });
  update();
})(typeof globalThis !== 'undefined' ? globalThis : this);
