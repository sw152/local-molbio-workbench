const test=require('node:test'),assert=require('node:assert/strict');
const map=require('../../src/localmolbio/static/sequence-map.js');
const f=(segments,label='feature')=>({label,type:'CDS',segments,renderable:true});
test('overlapping annotations never reuse a lane; adjacent disjoint annotations can',()=>{
 const features=[f([{start:0,end:40,strand:1}]),f([{start:20,end:60,strand:-1}]),f([{start:70,end:90,strand:1}])];
 const result=map.layout(features,100,false);
 assert.equal(result.lanes,2);
 for(const a of result.entries)for(const b of result.entries)if(a.id!==b.id&&a.lane===b.lane)assert.equal(map.overlaps(a.segments,b.segments,100,false),false);
});
test('origin-crossing multipart feature stays on one track and conflicts at both ends',()=>{
 const result=map.layout([f([{start:90,end:100,strand:-1},{start:0,end:15,strand:-1}]),f([{start:0,end:12}]),f([{start:93,end:99}])],100,true);
 const origin=result.entries.find(e=>e.id===0);
 assert.equal(origin.segments.length,2);
 assert.notEqual(origin.lane,result.entries.find(e=>e.id===1).lane);
 assert.notEqual(origin.lane,result.entries.find(e=>e.id===2).lane);
});
test('unavailable coordinates stay flagged, never clamped or invented',()=>{
 const result=map.layout([f([{start:-5,end:20}]),f([{start:20,end:20}]),f([{start:90,end:110}]),f([{start:0,end:10},{start:100,end:120}])],100,false);
 assert.equal(result.invalid,4);assert.equal(result.entries.length,1);assert.equal(result.entries[0].segments.length,1);
});
test('more than 60 nested annotations remain distinct',()=>{
 const result=map.layout(Array.from({length:70},(_,i)=>f([{start:i,end:1000-i,strand:i%2?1:-1}])),1000,true);
 assert.equal(result.entries.length,70);assert.equal(result.lanes,70);
});
test('arrow geometry respects strand and full-circle features remain finite',()=>{
 const forward=map.ribbon({start:20,end:40,strand:1},150,100,300);
 const reverse=map.ribbon({start:20,end:40,strand:-1},150,100,300);
 assert.notEqual(forward,reverse);
 for(const strand of [1,-1,null])assert.doesNotMatch(map.ribbon({start:0,end:100,strand},150,100,300),/NaN|Infinity/);
 assert.match(map.ribbon({start:0,end:100,strand:null},150,100,300),/circle/);
});
test('linear and unknown topology do not imply a circular molecule; labels are escaped',()=>{
 const features=[f([{start:0,end:100,strand:1}],'<script>unsafe</script>')];
 for(const topology of ['linear','unknown']) {
 const output=map.svgMarkup({features,length_bp:100,topology},map.layout(features,100,false));
 assert.doesNotMatch(output.markup,/<circle|<script>/);assert.match(output.markup,/<polygon/);assert.match(output.markup,/&lt;script&gt;/);
 }
});
test('ticks are ascending, bounded and distinct for small and large sequences',()=>{
 for(const length of [1,3,100,1001,10000]) {
 const ticks=map.ticks(length);assert.equal(ticks[0],0);assert(ticks.every(t=>t>=0&&t<length));assert.equal(new Set(ticks).size,ticks.length);
 }
});
