/* Direct Arrival: finite synthetic graph dynamics.
 * Two-mode rule: Engineering research exhibit, mechanism section 3.
 * Topology, gate order and angles are illustrative model inputs.
 * No physical post-isolation coupling is inferred. */
(function (root) {
  'use strict';
  const N = 6;
  function seed() {
    return Array.from({length: 2 * N}, (_, i) => [Math.sqrt(i % N === 0 ? .6 : .08), 0]);
  }
  function copy(state) { return state.map(z => z.slice()); }
  function account(state, start = 0, count = state.length) {
    return state.slice(start, start + count).reduce((sum, z) => sum + z[0] ** 2 + z[1] ** 2, 0);
  }
  function phase(state, node, angle) {
    const [r, i] = state[node], c = Math.cos(angle), s = Math.sin(angle);
    state[node] = [c * r - s * i, s * r + c * i];
  }
  // The same rule applies to EVERY edge, whether internal or between ledgers.
  // Old values are read before either endpoint is written.
  function exchange(state, [u, v, angle]) {
    const [ar, ai] = state[u], [br, bi] = state[v];
    const c = Math.cos(angle), s = Math.sin(angle);
    state[u] = [c * ar + s * bi, c * ai - s * br];
    state[v] = [c * br + s * ai, c * bi - s * ar];
  }
  function schedule(coupled) {
    const local = [];
    for (let group = 0; group < 2; group++) {
      for (let node = 0; node < N; node++) {
        local.push([group * N + node, group * N + (node + 1) % N, .035]);
      }
    }
    return [...local, [0, N, coupled ? .085 : 0], ...local.slice().reverse()];
  }
  function step(state, edges, reverse = false) {
    const sequence = reverse ? edges.slice().reverse().map(([u,v,t]) => [u,v,-t]) : edges;
    sequence.forEach(edge => exchange(state, edge));
  }
  function distance(a,b) { return Math.max(...a.flatMap((z,i) => z.map((x,k) => Math.abs(x-b[i][k])))); }
  root.DirectArrivalModel = {N, seed, copy, account, phase, exchange, schedule, step, distance};
  if (typeof module !== 'undefined' && module.exports) module.exports = root.DirectArrivalModel;
})(globalThis);
