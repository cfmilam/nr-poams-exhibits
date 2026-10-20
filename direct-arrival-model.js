/* View/controller for the separate stored graph states; no endpoint mirroring. */
(function(){
  'use strict';
  const M=globalThis.DirectArrivalModel, root=document.getElementById('portal-lab');
  const el=id=>root.querySelector('#pl-'+id), ns='http://www.w3.org/2000/svg';
  if(!root)return;
  if(!M){el('status').textContent='The interactive model could not load. Reload the page; the instructions and equations remain available.';return;}
  root.querySelectorAll('button,input').forEach(control=>{control.disabled=false;});
  let actual=M.seed(), control=M.seed(), history=[], tick=0, timer=null, trace=[{x:0,y:0}], message='Start with “Turn A +¼ cycle”, then “Evolve 10 steps”. B is not edited by the pulse.';
  let records=[];
  const signed=x=>(x>=0?'+':'')+x.toFixed(6);
  function svgNode(tag,attrs,text){const e=document.createElementNS(ns,tag);Object.entries(attrs||{}).forEach(([k,v])=>e.setAttribute(k,String(v)));if(text!==undefined)e.textContent=text;return e;}
  function graph(target,offset,color){
    target.replaceChildren();const points=Array.from({length:M.N},(_,i)=>{const a=-Math.PI/2+i*2*Math.PI/M.N;return [110+68*Math.cos(a),91+63*Math.sin(a)];});
    points.forEach(([x,y],i)=>{const [x2,y2]=points[(i+1)%M.N];target.append(svgNode('line',{x1:x,y1:y,x2,y2,stroke:'currentColor','stroke-opacity':.22,'stroke-width':2}));});
    points.forEach(([x,y],i)=>{const [r,im]=actual[offset+i],weight=r*r+im*im,phase=Math.atan2(im,r),radius=7+12*Math.sqrt(weight);
      target.append(svgNode('circle',{cx:x,cy:y,r:radius,fill:color,'fill-opacity':.18,stroke:color,'stroke-width':i===0?3:1.5}));
      target.append(svgNode('line',{x1:x,y1:y,x2:x+(radius-2)*Math.cos(phase),y2:y-(radius-2)*Math.sin(phase),stroke:color,'stroke-width':2}));
      target.append(svgNode('text',{x,y:y+radius+13,'text-anchor':'middle',fill:'currentColor','font-size':10},(offset?'B':'A')+i));
    });
  }
  function plot(){
    const target=el('trace');target.replaceChildren();const w=Math.max(280,target.getBoundingClientRect().width||640),h=150,left=42,right=12,top=16,bottom=30;
    target.setAttribute('viewBox','0 0 '+w+' '+h);
    const maxX=Math.max(10,...trace.map(p=>p.x)),maxY=Math.max(.1,...trace.map(p=>Math.abs(p.y)))*1.15;
    const xp=x=>left+x/maxX*(w-left-right),yp=y=>top+(1-y/maxY)/2*(h-top-bottom);
    target.append(svgNode('line',{x1:left,y1:yp(0),x2:w-right,y2:yp(0),stroke:'currentColor','stroke-opacity':.25}));
    target.append(svgNode('line',{x1:left,y1:top,x2:left,y2:h-bottom,stroke:'currentColor','stroke-opacity':.25}));
    [-maxY,0,maxY].forEach(y=>target.append(svgNode('text',{x:left-5,y:yp(y)+4,'text-anchor':'end',fill:'currentColor','font-size':11},y.toFixed(2))));
    [0,maxX].forEach(x=>target.append(svgNode('text',{x:xp(x),y:h-10,'text-anchor':x?'end':'start',fill:'currentColor','font-size':11},String(x))));
    target.append(svgNode('polyline',{points:trace.map(p=>xp(p.x)+','+yp(p.y)).join(' '),fill:'none',stroke:'var(--violet)','stroke-width':2.5}));
  }
  function response(){return M.account(actual,M.N,M.N)-M.account(control,M.N,M.N);}
  function render(){
    graph(el('a'),0,'var(--blue)');graph(el('b'),M.N,'var(--violet)');
    el('qa').textContent=M.account(actual,0,M.N).toFixed(6);el('qb').textContent=M.account(actual,M.N,M.N).toFixed(6);
    el('response').textContent=signed(response());el('total').textContent=M.account(actual).toFixed(6);el('error').textContent=Math.abs(M.account(actual)-2).toExponential(1);
    el('tick').textContent='Step '+tick;el('status').textContent=message;el('run').textContent=timer?'Pause':'Run';el('undo').disabled=history.length===0;
    el('bridge').textContent=el('link').checked?'A0 ↔ B0 · explicit coupling assumption':'No exchange edge · two independent graphs';el('bridge').classList.toggle('off',!el('link').checked);
    el('log').textContent=records.slice(-3).join('\n')||'Waiting for an operation.';plot();
  }
  function pause(){if(timer){clearInterval(timer);timer=null;} }
  function pulse(angle){
    pause();const before=M.copy(actual);M.phase(actual,0,angle);history.push({type:'phase',angle,tick});
    const bUnchanged=M.distance(before.slice(M.N),actual.slice(M.N))===0;
    message='A0 phase changed. B’s stored state is '+(bUnchanged?'unchanged':'unexpectedly changed')+' immediately. Evolve to test the interaction.';
    records.push('A-only pulse '+(angle>0?'+':'−')+'¼ cycle; immediate B change = 0.');render();
  }
  function evolve(count=1){
    for(let i=0;i<count;i++){
      const edges=M.schedule(el('link').checked), beforeA=M.account(actual,0,M.N),beforeB=M.account(actual,M.N,M.N);
      M.step(actual,edges);M.step(control,edges);history.push({type:'step',edges,tick});tick++;
      const da=M.account(actual,0,M.N)-beforeA,db=M.account(actual,M.N,M.N)-beforeB;
      records.push('Step '+tick+': ΔA '+signed(da)+'  ΔB '+signed(db)+'  sum '+(da+db).toExponential(1));trace.push({x:tick,y:response()});
    }
    message=el('link').checked?'The same edge law evolves both graphs. B’s response versus the unpulsed run is '+signed(response())+'.':'The edge is off. No new A→B influence is admitted; any earlier transferred state remains. Reset for a clean isolation test.';
    render();
  }
  el('plus').addEventListener('click',()=>pulse(Math.PI/2));el('minus').addEventListener('click',()=>pulse(-Math.PI/2));
  el('step').addEventListener('click',()=>{pause();evolve(10);});
  el('run').addEventListener('click',()=>{if(timer){pause();render();return;}let steps=0;timer=setInterval(()=>{evolve();if(++steps>=200){pause();render();}},100);render();});
  el('link').addEventListener('change',()=>{pause();message='Exchange setting changed for future steps; stored states have not been erased.';render();});
  el('reset').addEventListener('click',()=>{pause();actual=M.seed();control=M.seed();history=[];records=[];tick=0;trace=[{x:0,y:0}];message='Reset complete. Turn A, then evolve. The exchange checkbox keeps your chosen setting.';render();});
  el('undo').addEventListener('click',()=>{
    pause();const op=history.pop();if(!op)return;
    if(op.type==='step'){M.step(actual,op.edges,true);M.step(control,op.edges,true);tick=op.tick;trace.pop();}else{M.phase(actual,0,-op.angle);}
    records.push('Undid '+op.type+' using its inverse operation, not a saved state.');message='Last operation reversed. The inverse acts on the actual stored state.';render();
  });
  if(typeof ResizeObserver!=='undefined'){const resizeObserver=new ResizeObserver(()=>plot());resizeObserver.observe(el('trace'));}
  else{window.addEventListener('resize',plot);}
  document.addEventListener('visibilitychange',()=>{if(document.hidden){pause();render();}});
  render();
})();
