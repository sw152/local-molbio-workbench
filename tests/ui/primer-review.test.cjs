const {test}=require('node:test'),assert=require('node:assert/strict');
const {pairsFrom,diagram}=require('../../src/localmolbio/static/primer-review.js');
const make=(direction,job='run')=>({direction,name:direction,sequence_text:'ACGT',binding_start:direction==='forward'?10:90,binding_end:direction==='forward'?30:110,metrics:{analysis_job_id:job,pair_index:1,product_size:100,target:{start:40,end:80}}});
test('pair identity includes analysis job and never joins unrelated candidates',()=>{
 assert.equal(pairsFrom([make('forward','a'),make('reverse','b')]).length,0);
 assert.equal(pairsFrom([make('forward','a'),make('reverse','a'),make('forward','b'),make('reverse','b')]).length,2);
});
test('invalid, duplicate and legacy candidates cannot create plausible diagrams',()=>{
 const f=make('forward'),r=make('reverse');assert.equal(pairsFrom([f,r,f]).length,0);
 assert.equal(pairsFrom([f,{...r,binding_end:111}]).length,0);
 assert.equal(pairsFrom([{...f,metrics:{}},r]).length,0);
 assert.equal(pairsFrom([f,{...r,metrics:{...r.metrics,target:{start:5,end:80}}}]).length,0);
});
test('diagram expresses true window coordinates, target and opposing arrows',()=>{
 const svg=diagram(pairsFrom([make('forward'),make('reverse')])[0]);
 assert.match(svg,/100 bp expected product/);assert.match(svg,/Target 41–80 bp/);
 assert.match(svg,/60,103 179,103 188,112/);assert.match(svg,/700,103 581,103 572,112/);
 assert.match(svg,/>110<\/text>/);assert.doesNotMatch(svg,/NaN|undefined/);
});

const {productReview}=require('../../src/localmolbio/static/primer-review.js');
test('missing or truncated evidence never claims no alternative products',()=>{
 const pair={left:{metrics:{}}};assert.match(productReview(pair),/not recorded/);
 pair.left.metrics.product_review={site_search_complete:false,observed_alternative_count:0,observed_product_count:0,product_size_min:100,product_size_max:800,products:[]};
 const html=productReview(pair);assert.match(html,/Incomplete site search/);assert.doesNotMatch(html,/No additional product/);
 pair.left.metrics.product_review.site_search_complete=true;
 assert.match(productReview(pair),/No additional product interval found in this range/);
});
