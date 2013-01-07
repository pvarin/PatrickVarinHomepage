function addGridItemEvents(){
	overlayItems = $('.grid-item-overlay');
	overlayItems.hide();
	
	hideOverlay = function (event){
		$(event.currentTarget).children('.grid-item-overlay').fadeOut('fase');
	};
	
	showOverlay = function (event){
		$(event.currentTarget).children('.grid-item-overlay').fadeIn('fast');
	};
	
	for (var i=0;i<overlayItems.length;i++){
		var parent = $(overlayItems[i]).parent()
		parent.hover(showOverlay,hideOverlay);
	}
}


$(document).ready(addGridItemEvents)
