function showOverlay(item){
	item.getElementsByClassName('grid-item-overlay-hidden')[0].setAttribute('class','grid-item-overlay-visible');
}

function hideOverlay(item){
	item.getElementsByClassName('grid-item-overlay-visible')[0].setAttribute('class','grid-item-overlay-hidden');
}

function addGridItemEvents(){
	overlayItems = document.getElementsByClassName('grid-item-overlay-hidden');
	console.log(overlayItems);
	for (var i=0; i<overlayItems.length; i++){
		overlayItems[i].parentNode.setAttribute('onmouseout','hideOverlay(this);');
		overlayItems[i].parentNode.setAttribute('onmouseover','showOverlay(this);');
	}
}
