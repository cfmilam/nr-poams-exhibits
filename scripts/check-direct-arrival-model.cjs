'use strict';
const assert=require('node:assert/strict');
const M=require('../direct-arrival-model-core.js');
const checks=[];
function check(name,fn){const detail=fn();checks.push({name,passed:true,...detail});}
check('Source-only pulse leaves every B component unchanged',()=>{
  const s=M.seed(),before=M.copy(s);M.phase(s,0,Math.PI/2);
  assert.deepEqual(s.slice(6),before.slice(6));assert(M.distance(s,before)>.1);
  return {maxBChange:M.distance(s.slice(6),before.slice(6))};
});
check('Fresh source intervention changes B only with exchange enabled',()=>{
  const a=M.seed(),ref=M.seed();M.phase(a,0,Math.PI/2);
  for(let i=0;i<10;i++){M.step(a,M.schedule(true));M.step(ref,M.schedule(true));}
  const response=M.account(a,6,6)-M.account(ref,6,6);assert(Math.abs(response)>.01);
  return {steps:10,receiverResponse:response,totalAccount:M.account(a)};
});
check('Disconnected graphs remain causally independent under many late source pulses',()=>{
  const a=M.seed(),ref=M.seed();
  for(let i=0;i<500;i++){M.phase(a,i%6,.37+i*.001);M.step(a,M.schedule(false));M.step(ref,M.schedule(false));}
  assert.deepEqual(a.slice(6),ref.slice(6));return {steps:500,maxBChange:0};
});
check('Arbitrary states preserve norm through long evolution',()=>{
  let maxError=0;
  for(let trial=0;trial<12;trial++){
    const a=Array.from({length:12},(_,i)=>[Math.sin(i*1.71+trial),Math.cos(i*.39-trial)]);
    const before=M.account(a);
    for(let i=0;i<1000;i++){M.step(a,M.schedule(i%3!==0));if(i%17===0)M.phase(a,0,-.8);}
    maxError=Math.max(maxError,Math.abs(M.account(a)-before));
  }
  assert(maxError<1e-9);return {trials:12,stepsPerTrial:1000,maxAbsoluteNormError:maxError};
});
check('Full complex state recovered by inverse sequence without snapshots',()=>{
  const a=M.seed(),initial=M.copy(a),ops=[];
  for(let i=0;i<300;i++){const edges=M.schedule(i%4!==0);M.step(a,edges);ops.push({edges});if(i%11===0){M.phase(a,0,.4);ops.push({phase:.4});}}
  ops.reverse().forEach(op=>op.edges?M.step(a,op.edges,true):M.phase(a,0,-op.phase));
  const error=M.distance(a,initial);assert(error<1e-11);return {maxComponentError:error};
});
check('Exchange rule independently matches single-edge analytic prediction',()=>{
  const phi=.73,theta=.29;const a=[[Math.cos(phi),Math.sin(phi)],[1,0]];
  M.exchange(a,[0,1,theta]);const expected=Math.sin(2*theta)*Math.sin(phi),actual=M.account(a,1,1)-1;
  assert(Math.abs(actual-expected)<1e-14);return {predictedReceiverChange:expected,computedReceiverChange:actual};
});
check('Identical graphs without intervention give zero reference response',()=>{
  const a=M.seed(),b=M.seed();assert.notEqual(a,b);assert.notEqual(a[0],b[0]);
  for(let i=0;i<100;i++){M.step(a,M.schedule(true));M.step(b,M.schedule(true));}
  assert.deepEqual(a,b);return {steps:100,referenceDifference:0};
});
console.log(JSON.stringify(checks,null,2));
