<!-- Global site tag (gtag.js) - Google Analytics -->
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  gtag('config', 'UA-111320548-1');


$(document).on('shiny:inputchanged', function(event) {
    if (event.name.substring(1,1)=='.') {
	gtag('event', 'search', { 'event_category': event.name, 'event_label': event.value } );
    }
});
