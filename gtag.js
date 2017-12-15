<!-- Global site tag (gtag.js) - Google Analytics -->
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  gtag('config', 'UA-111320548-1');


$(document).on('shiny:inputchanged', function(event) {
    gtag('event', 'input', { name: event.name, value: event.value } );
});
