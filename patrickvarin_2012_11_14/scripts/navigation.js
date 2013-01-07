/*Navigation Functions*/
currentPage = 'home'
function goToPage(page){
	//hide all
	$('#'+currentPage).fadeOut(function(){
	$('#'+page).fadeIn();});
	//show content
	currentPage = page;

	//remove navbar formatting
	$('#navbar li').not($('#'+page+'-button')).toggleClass('current-tab',false);
	//format navbar
	$('#'+page+'-button').toggleClass('current-tab',true);
}

function goToGrid(grid){
	//hide all
	document.getElementById('projects').setAttribute("class",'content-grid-wrapper-hidden');
	document.getElementById('research').setAttribute("class",'content-grid-wrapper-hidden');
	document.getElementById('art').setAttribute("class",'content-grid-wrapper-hidden');
	
	//remove menu formatting
	document.getElementById('projects-button').setAttribute("class","");
	document.getElementById('research-button').setAttribute("class","");
	document.getElementById('art-button').setAttribute("class","");
	
	//show content
	document.getElementById(grid).setAttribute("class","content-grid-wrapper-visible")
	
	//format navbar
	document.getElementById(grid+'-button').setAttribute('class','current-grid')
}

$(document).ready(function(){
	$('.content-wrapper').not($('#home')).hide();
});
