function goToPage(page){
	//hide all
	document.getElementById('home').setAttribute("class","content-wrapper-hidden");
	document.getElementById('about').setAttribute("class","content-wrapper-hidden");
	document.getElementById('grid').setAttribute("class","content-wrapper-hidden");
	document.getElementById('resume').setAttribute("class","content-wrapper-hidden");
	document.getElementById('contact').setAttribute("class","content-wrapper-hidden");

	//remove navbar formatting
	document.getElementById('home-button').setAttribute("class","");
	document.getElementById('about-button').setAttribute("class","");
	document.getElementById('grid-button').setAttribute("class","");
	document.getElementById('resume-button').setAttribute("class","");
	document.getElementById('contact-button').setAttribute("class","");
	
	//show content
	document.getElementById(page).setAttribute("class","content-wrapper-visible");
	
	//format navbar
	document.getElementById(page+'-button').setAttribute("class","current-tab");
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
