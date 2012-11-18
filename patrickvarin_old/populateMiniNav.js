// JavaScript Document
navList = document.getElementById('mininav');
entryList = document.getElementsByClassName('content-entry');
for (i=0;i<entryList.length;i++){
	listItem = document.createElement('li');
	linkItem = document.createElement('a');
	linkItem.innerHTML = entryList[i].id;
	linkItem.href = '#'+entryList[i].id;
	listItem.appendChild(linkItem);
	navList.appendChild(listItem);
}