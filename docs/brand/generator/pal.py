import colorsys
def _lin(c):
    c=c/255; return c/12.92 if c<=0.03928 else ((c+0.055)/1.055)**2.4
def lum(h):
    h=h.lstrip('#'); r,g,b=(int(h[i:i+2],16) for i in (0,2,4)); return 0.2126*_lin(r)+0.7152*_lin(g)+0.0722*_lin(b)
def cr(a,b):
    la,lb=sorted([lum(a),lum(b)],reverse=True); return (la+0.05)/(lb+0.05)
def hexs(r,g,b): return '#%02x%02x%02x'%tuple(max(0,min(255,round(v*255))) for v in (r,g,b))
def shift(h,dl):
    h=h.lstrip('#'); r,g,b=(int(h[i:i+2],16)/255 for i in (0,2,4))
    H,L,S=colorsys.rgb_to_hls(r,g,b); L=max(0,min(1,L+dl)); return hexs(*colorsys.hls_to_rgb(H,L,S))
def fit(h,ground,target,direction):
    while cr(h,ground)<target: h=shift(h,0.01*direction)
    return h
WASHI=dict(ground="#f3ede1",alt="#ebe3d3",raised="#f7f2e8",line="#ddd4c3",ink="#1d1b18",muted="#5f5a52",mesh="#6f6a62")
NIGHT=dict(ground="#0e1116",alt="#0a0d11",raised="#141820",line="#222831",ink="#e8ecf0",muted="#a9b1bb",mesh="#8b949e")
RAW={
 "tokusa":("Tokusa","木賊","horsetail green, the colour of moss between garden stones","#2e6b55","#5fd3a3"),
 "yamabuki":("Yamabuki","山吹","kerria gold, autumn leaves on a bamboo grove","#9a6200","#ffc23a"),
 "ruri":("Ruri","瑠璃","lapis, the deep blue of Hiroshige's water","#1f4fa6","#7fa9ff"),
 "benibana":("Benibana","紅花","safflower crimson, the red of a monk's lacquer bowl","#b0243f","#ff7f96"),
 "kikyo":("Kikyō","桔梗","bellflower violet, moonlight through paper","#5a3ea3","#b8a2ff"),
 "hisui":("Neon Hisui","翡翠","electric jade, an abacus bead lit from inside","#007a68","#2dffd5"),
}
PALS={}
for k,(n,kj,desc,a,b) in RAW.items():
    acc=fit(a,WASHI["ground"],4.8,-1); link=fit(shift(acc,-0.04),WASHI["ground"],5.5,-1)
    nacc=fit(b,NIGHT["ground"],7.0,+1); nlink=fit(shift(nacc,0.04),NIGHT["ground"],8.0,+1)
    PALS[k]=dict(key=k,name=n,kanji=kj,desc=desc,
        washi=dict(WASHI,acc=acc,link=link,cr=round(cr(acc,WASHI["ground"]),1),crl=round(cr(link,WASHI["ground"]),1)),
        night=dict(NIGHT,acc=nacc,link=nlink,cr=round(cr(nacc,NIGHT["ground"]),1),crl=round(cr(nlink,NIGHT["ground"]),1)))
if __name__=="__main__":
    for k,p in PALS.items(): print(k,p["washi"]["acc"],p["washi"]["cr"],p["night"]["acc"],p["night"]["cr"])
