document.addEventListener('DOMContentLoaded',()=>{
  const button=document.querySelector('.menu-button');
  const nav=document.querySelector('.sidebar');
  if(button&&nav){button.addEventListener('click',()=>{const open=nav.classList.toggle('open');button.setAttribute('aria-expanded',String(open));});}
  const search=document.querySelector('#appendix-search');
  const cards=[...document.querySelectorAll('.appendix-card')];
  const empty=document.querySelector('#no-results');
  if(search){search.addEventListener('input',()=>{const q=search.value.trim().toLowerCase();let shown=0;cards.forEach(card=>{const show=!q||card.dataset.search.includes(q);card.hidden=!show;if(show)shown++;});if(empty)empty.hidden=shown!==0;});}
});
