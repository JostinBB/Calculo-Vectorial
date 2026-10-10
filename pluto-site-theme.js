/* Sin dependencias: usa el mismo tema para el índice y los notebooks. */
(function () {
  'use strict';
  var root = document.documentElement;
  var toggle = document.querySelector('[data-cv-theme-toggle]');
  function updateButton() {
    if (!toggle) return;
    var dark = root.dataset.theme !== 'light';
    toggle.setAttribute('aria-label', dark ? 'Activar modo claro' : 'Activar modo oscuro');
    toggle.setAttribute('title', dark ? 'Activar modo claro' : 'Activar modo oscuro');
    toggle.querySelector('[data-cv-icon]').innerHTML = dark
      ? '<circle cx="12" cy="12" r="4"/><path d="M12 2v2m0 16v2M4.93 4.93l1.42 1.42m11.3 11.3 1.42 1.42M2 12h2m16 0h2M4.93 19.07l1.42-1.42m11.3-11.3 1.42-1.42"/>'
      : '<path d="M20.7 14.1A8.5 8.5 0 0 1 9.9 3.3 8.6 8.6 0 1 0 20.7 14.1Z"/>';
  }
  if (toggle) {
    updateButton();
    toggle.addEventListener('click',function () {
      root.dataset.theme = root.dataset.theme === 'dark' ? 'light' : 'dark';
      try { localStorage.setItem('taller-theme', root.dataset.theme); } catch (_) {}
      updateButton();
    });
  }
})();
