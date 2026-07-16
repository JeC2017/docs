window.ftsMathJaxReady=(async()=>{
  const inject=text=>{const script=document.createElement('script');script.textContent=text;document.head.appendChild(script);};
  const fallback=()=>new Promise((resolve,reject)=>{const script=document.createElement('script');script.src='https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js';script.onload=resolve;script.onerror=reject;document.head.appendChild(script);});
  try{
    if(typeof DecompressionStream==='undefined')throw new Error('gzip streaming unavailable');
    const response=await fetch('assets/vendor/mathjax-3.2.2-tex-svg.js.gz');
    if(!response.ok||!response.body)throw new Error('local MathJax asset unavailable');
    const stream=response.body.pipeThrough(new DecompressionStream('gzip'));
    inject(await new Response(stream).text());
  }catch(error){await fallback();}
  if(window.MathJax&&window.MathJax.startup&&window.MathJax.startup.promise){await window.MathJax.startup.promise;}
  return true;
})();
