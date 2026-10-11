/* Fixed corpus reading theme. Run before legacy page/canvas initialization. */
'use strict';
document.documentElement.dataset.theme = 'dark';
try {
  localStorage.setItem('poams-theme', 'dark');
  localStorage.setItem('poamsTheme', 'dark');
  localStorage.setItem('poams-analysis-theme', 'dark');
  localStorage.setItem('navigation-engine-theme', 'dark');
} catch (_) { /* Reading works when browser storage is unavailable. */ }
