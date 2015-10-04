/// Tools for oligonucleotide design
namespace Amyris
open System

/// Core primer Tm calculation
module tmcore =
    [<Measure>] type M // Molar
    [<Measure>] type uM // micromolar
    [<Measure>] type mM // millimolar
    [<Measure>] type nM // nanomolar
    [<Measure>] type K // kelvin
    [<Measure>] type C // Celsius
    [<Measure>] type Cal // Calories
    let inline K2C (k:float<K>) = (k*1.0<C/K>) - 273.5<C>
    let inline C2K (c:float<C>) = (c*1.0<K/C>) + 273.5<K>
    let inline mM2M (c:float<mM>) = c*0.001<M/mM>
    let inline uM2M (c:float<uM>) = c/1000000.0<uM/M>
    let inline nM2M (c:float<nM>) = c/1000000000.0<nM/M> // Table 2: http://pubs.acs.org/doi/pdf/10.1021/bi702363u
    let private a = 3.92E-5<1/K>
    let private b = -9.11E-6<1/K>
    let private c = 6.26E-5<1/K>
    let private d = 1.42E-5<1/K>
    let private e = -4.82E-4<1/K>
    let private f = 5.25E-4<1/K>
    let private g = 8.31E-5<1/K>

    /// fGC - proportion of GC bases in full primer array
    let gc (primer : char array) = (Array.fold (fun gc c -> match c with | 'G' | 'g' | 'C' | 'c' -> gc+1 | _ -> gc) 0 primer |> float) / (float primer.Length)
    
    /// fGC - proportion of GC bases in first N bases of primer
    let gcN (primer : char array) N = 
        let rec count i gc =
            if i = N then  (float gc)/(float N)
            else match primer.[i] with | 'G' | 'g' | 'C' | 'c' ->  count (i+1) (gc+1) | _ -> count (i+1) gc

        count 0 0

    /// Take natural log of a molar constant producing a unitless value
    let inline logMolar (x:float<M>) = log ( x / 1.0<M>)

    //let targetLength = 20
    // Previous distance based penalty for not ending on a G or C
    //let ATPenalty = 5

    ///<summary>Parameters used to design primers. These include penalties used to 
    ///rank primers using the various criteria and settings for the calculation of Tm. 
    ///</summary>
    ///<typeparam name="lengthPenalty"> Penalty parameter for difference from targeted primer length. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="tmPenalty"> Penalty parameter for difference from targeted Tm. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="positionPenalty"> Penalty parameter for difference from targeted position. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="polyLengthThreshold"> Upper admissible length for polynucleotide. </typeparam> 
    ///<typeparam name="polyPenalty"> Penalty for going beyong the polynucleotide length limit. </typeparam> 
    ///<typeparam name="threePrimeUnstablePenalty"> Penalty for having an unstable 3' end. </typeparam> 
    ///<typeparam name="ATPenalty"> Penalty for the primer ending on a A or T. </typeparam> 
    ///<typeparam name="maxLength"> Maximum primer length. </typeparam> 
    ///<typeparam name="minLength"> Minimum primer length. </typeparam> 
    ///<typeparam name="targetLength"> Target primer length. </typeparam> 
    ///<typeparam name="monovalentConc"> Monovalent cation concentration (Na+, K+) in mol/L. </typeparam> 
    ///<typeparam name="divalentConc"> Divalent cation concentration (Mg2+, Ca2+) in mol/L. </typeparam> 
    ///<typeparam name="primerConc"> Primer concentration in mol/L. </typeparam> 
    ///<typeparam name="templateConc"> Template concentration in mol/L. </typeparam> 
    ///<typeparam name="dNTPConc"> dNTP concentration in mol/L. </typeparam> 
    type PrimerParams = {   lengthPenalty : float ; 
                            tmPenalty : float ; 
                            tmMaxDifference : float<C> ;
                            positionPenalty : float ; 
                            polyLengthThreshold : int ; 
                            polyPenalty : float ; 
                            threePrimeUnstablePenalty : float ; 
                            ATPenalty :float<C> ;
                            maxLength : int ; 
                            minLength : int ; 
                            targetLength : int ;
                            monovalentConc : float<M> ; 
                            divalentConc : float<M> ; 
                            primerConc : float<M> ; 
                            templateConc : float<M> ; 
                            dNTPConc : float<M>}
    
    ///<summary>Parameters used to design primers. These include penalties used to 
    ///rank primers using the various criteria and settings for the calculation of Tm. 
    ///</summary>
    ///<typeparam name="lengthPenalty"> Penalty parameter for difference from targeted primer length. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="tmPenalty"> Penalty parameter for difference from targeted Tm. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="positionPenalty"> Penalty parameter for difference from targeted position. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="polyLengthThreshold"> Upper admissible length for polynucleotide. </typeparam> 
    ///<typeparam name="polyPenalty"> Penalty for going beyong the polynucleotide length limit. </typeparam> 
    ///<typeparam name="threePrimeUnstablePenalty"> Penalty for having an unstable 3' end. </typeparam> 
    ///<typeparam name="ATPenalty"> Penalty for the primer ending on a A or T. </typeparam> 
    ///<typeparam name="maxLength"> Maximum primer length. </typeparam> 
    ///<typeparam name="minLength"> Minimum primer length. </typeparam> 
    ///<typeparam name="targetLength"> Target primer length. </typeparam> 
    ///<typeparam name="monovalentConc"> Monovalent cation concentration (Na+, K+) in mol/L. </typeparam> 
    ///<typeparam name="divalentConc"> Divalent cation concentration (Mg2+, Ca2+) in mol/L. </typeparam> 
    ///<typeparam name="primerConc"> Primer concentration in mol/L. </typeparam> 
    ///<typeparam name="templateConc"> Template concentration in mol/L. </typeparam> 
    ///<typeparam name="dNTPConc"> dNTP concentration in mol/L. </typeparam> 
    let defaultParams = {   tmPenalty = 1.0; 
                            tmMaxDifference = 5.0<C> ;
                            positionPenalty = 5.0 ; 
                            lengthPenalty = 3.0 ; 
                            polyLengthThreshold = 4; 
                            polyPenalty = 10. ; 
                            threePrimeUnstablePenalty = 5.0 ; 
                            ATPenalty=3.0<C> ; 
                            targetLength=20 ;
                            maxLength = 60 ;
                            minLength = 20 ; 
                            monovalentConc = mM2M 50.0<mM>;
                            primerConc = uM2M 0.25<uM> ; 
                            divalentConc = mM2M 1.5<mM> ; 
                            templateConc = uM2M 0.01<uM> ; 
                            dNTPConc = uM2M 0.0<uM> ;}

    
    /// Final oligo design result with the oligo, its melting temp and offset into the template
    type Oligo = { tag : string ; oligo : char array ; temp : float<C> ; offset : int }

    type GeneInfo = { systematic : string ; name : string ; sequence : char [] }

    type ANCHOR = LEFT | RIGHT | CENTERLEFT | CENTERRIGHT
    type STRAND = TOP | BOTTOM

    /// Inputs to the primer designer. 
    /// temp: template sequence
    /// align: how to align the primer with respect to the template
    /// strand: on what strand to search the primer
    /// offset: 
    /// sequencePenalties: a float [] option of the same size as temp, giving nucleotide specific penalties. If None, an array of 0. will be used.
    ///    These penalties are only used in the CENTERLEFT and CENTERRIGHT ANCHOR types.
    type OligoTask = { tag : string ; temp : char [] ; align :ANCHOR ; strand : STRAND ; offset : int ; targetTemp : float<C> ; sequencePenalties: float [] option}
    
    type OligoDesign = { tag : string ; oligo : string ; temp : float<C> ; note : string ; pos : int}

    type PrimerHit = { chr : int ; sysName : string ; tag : string ; pos : int ; fwd : bool}
    type PCRHit = { p1 : PrimerHit ; p2 : PrimerHit }
   
    let base2GC b =
        match b with
            | 'G' | 'C' | 'g' | 'c' -> 1
            | 'A' | 'a' | 'T' | 't' -> 0
            | _ -> failwith (sprintf "bad base %c in base2GC" b)

    type HS = { h : float<Cal/M> ; s: float<Cal/K/M>}
    /// Santa lucia PNAS Feb 17 1998 vol 95 no 3. 1460-1465
    let HSpair (a:char) (b:char) : HS =
            match a,b with
                | 'A','A' | 'T','T'     -> {h= -7.9<Cal/M>; s= -22.2<Cal/M/K> }
                | 'A','G' | 'C','T'     -> {h= -7.8<Cal/M>; s= -21.0<Cal/M/K>}
                | 'A','T'               -> {h= -7.2<Cal/M>; s= -20.4<Cal/M/K> }
                | 'A','C'| 'G','T'      -> {h= -8.4<Cal/M>; s= -22.4<Cal/M/K>}
                | 'G','A' | 'T','C'     -> {h= -8.2<Cal/M>; s= -22.2<Cal/M/K> }
                | 'G','G' | 'C','C'     -> {h= -8.0<Cal/M>; s= -19.9<Cal/M/K> }
                | 'G','C'               -> {h= -9.8<Cal/M>; s= -24.4<Cal/M/K>}
                | 'T','A'               -> {h= -7.2<Cal/M>; s= -21.3<Cal/M/K> }
                | 'T','G' | 'C','A'     -> {h= -8.5<Cal/M>; s= -22.7<Cal/M/K>}
                | 'C','G'               -> {h= -10.6<Cal/M>; s= -27.2<Cal/M/K>}
                | _ -> failwith "bad nucpair"

    let deltaG (a:char) (b:char) =
        match a,b with
            | 'A','A' | 'T','T'     -> -1.0
            | 'A','T'               -> -0.88
            | 'T','A'               -> -0.58
            | 'C','A' | 'T','G'     -> -1.45
            | 'G','T' | 'A','C'     -> -1.44
            | 'C','T' | 'A','G'     -> -1.28
            | 'G','A' | 'T','C'     -> -1.30
            | 'C','G'               -> -2.17
            | 'G','C'               -> -2.24
            | 'G','G' | 'C','C'     -> -1.84
            | _ -> sprintf "bad nucpair /%c/%c/" a b |> failwith    

    let deltaG2 (a:char) (b:char) =
        match a,b with
            | 'A','A' | 'T','T'     -> -1.9
            | 'A','T'               -> -1.5
            | 'T','A'               -> -0.9
            | 'C','A' | 'T','G'     -> -1.9
            | 'G','T' | 'A','C'     -> -1.3
            | 'C','T' | 'A','G'     -> -1.6
            | 'G','A' | 'T','C'     -> -1.6
            | 'C','G'               -> -3.6
            | 'G','C'               -> -3.1
            | 'G','G' | 'C','C'     -> -3.1
            | _ -> sprintf "bad nucpair /%c/%c/" a b |> failwith   
    let initG2 (c:char) =
            match c with
                | 'G' | 'C' -> -1.96 // 5.0
                | 'A' | 'T' -> -1.96 - 0.05 // 6.0
                | _  -> sprintf "bad nuc %c in initG2" c |> failwith
                     
    /// Test if an oligo is self complementary
    let selfComp (oligo:char array) length =
        let rec checkOne l r =
            l>r ||
                match oligo.[l],oligo.[r] with
                    | 'G','C' | 'C','G' | 'A','T' | 'T','A' -> checkOne (l+1) (r-1)
                    | _ -> false
        if length%2=1 then false else checkOne 0 (length-1)

    /// NN implementation with 1M Na+ concentration
    /// http://www.pnas.org/content/95/4/1460.full
    /// http://en.wikipedia.org/wiki/DNA_melting
    let wikiTemp1000 (primerA:float<M>) (primerB:float<M>) (o : char array) N =
        // Nucleic Acids Research, 1996, Vol. 24, No. 22 4501–4505
        let endBase1 (a:char) = match a with 'A' | 'T' -> {h=2.3<Cal/M> ; s=4.1<Cal/M/K>}  | 'C' | 'G' -> {h=0.1<Cal/M>;s= -2.8<Cal/M/K>} | _ -> failwith (sprintf "bad start base '%c'" a)

        let rec score i hs =
            if i = N then 
                let hs' = endBase1 (o.[i-1])
                {h=hs.h + hs'.h ; s = hs.s + hs'.s}
            else
                let hs' = HSpair (o.[i-1]) (o.[i]) 
                score (i+1) {h= hs.h+hs'.h ; s= hs.s + hs'.s } // +(h+h') (s+s')
        let hSF = endBase1 (o.[0]) // 0.0,0.0 
        let hs = 
            match N with
                | 1 -> hSF // hS,hF
                | 2 ->
                    let hSF2 = endBase1 (o.[1])  // Start/finish
                    {h=hSF.h+hSF2.h ;  s= hSF.s+hSF2.s}
                | _ -> score 1 hSF
        let R = 1.98722<Cal/M/K> // ideal gas constant

        // Santalucia
        // Proc.Natl.Acad. Sci.USA Vol. 95, pp. 1460–1465, February 1998 Biochemistry
        // eq (3)
        // Self complementary molecules get Ct = primerConcentration
        // Non self complementary case gets primerConc/4 if strands are in equal concentration
        // or (Ca-Cb/2) if different concentrations were Ca and Cb strands respectively and Ca > Cb
        let comparable (a:float<M>) (b:float<M>) = abs(a-b)/(a+b) < 0.1
        let symmetrical = selfComp o N 
        let Ct = if symmetrical then primerA else
                    if comparable primerA primerB then primerA / 4.0  else ( (max primerA primerB)-(min primerA primerB)/2.0)
               
        let hs' = { h = hs.h * 1000.0 ; s = hs.s } // This magic is due to deltaH units being in 100Cal/mol and deltaS being in 0.1Cal/Mol 
        hs'.h/(
                ( hs'.s+(if symmetrical then -1.4<Cal/K/M> else 0.0<Cal/K/M>)+(R*logMolar(Ct) ) )
             ) 
                 

    /// Calc temp given a GC and oligo length
    let gcN2Temp gc N = 
        if N < 14 then (4*gc) + 2 * (N-gc)*2
        else (64.9+41.0*(float(gc)-16.4)/float(N) ) |> int

    /// Monovalent salt correction
    let monoSaltCorrect (mon:float<M>) (Ca:float<M>) (Cb:float<M>) (o : char array) N =
        let tm1000 = wikiTemp1000 Ca Cb o N
        let fgc = gcN o N
        //http://www.owczarzy.net/Tm_for_duplexes-IDT_Tech.pdf
        let correction = (4.29*fgc-3.95)*1E-5<1/K>*logMolar(mon)+9.40E-6<1/K>*logMolar(mon)*logMolar(mon)
        let tmMonRecip = (1.0/tm1000) + correction
        1.0/tmMonRecip |> K2C

    // eq 16 table 2
    let eq16Table2 (div:float<M>) (primerA:float<M>) (primerB:float<M>) (dNTP:float<M>) (primer: char array) (N:int) =
        let mgConc = max 0.0<M> (div-dNTP) // NTPs bing divalent ions, so remove that from further consideration
        let l = logMolar mgConc
        let tm1000 = wikiTemp1000 primerA primerB primer N
    
        // fGC is the fraction of residues that are G or C
        let fGC = gcN primer N

        let mgRecip = 1.0<1> /tm1000 + a + b*l + fGC * (c+d*l) + (1.0/(2.0 * float (N-1)) * (e+f*l+g*l*l) )
        1.0/mgRecip |> K2C

    // eq 16 table 2
    let eq16Table2Mod (div:float<M>) (mon:float<M>) (Ca:float<M>) (Cb:float<M>) (dNTP:float<M>) (primer: char array) (N:int) =
        let mgConc = max 0.0<M> (div-dNTP) // NTPs bind divalent ions, so remove that from further consideration
        let l = logMolar mgConc
        let lm = logMolar mon
        let tm1000 = wikiTemp1000 Ca Cb primer N
    
        let modA  = 3.92E-5<1/K> * (0.843 - 0.352 * sqrt(mon / 1.0<M> (* UNITHACK *) ) * lm)
        let modD  = 1.42E-5<1/K> *(1.279-4.03E-3*lm-8.03E-3*lm*lm)
        let modG  = 8.31E-5<1/K> *(0.486-0.258*lm+5.25E-3*lm*lm*lm)

        // fGC is the fraction of residues that are G or C
        let fGC = gcN primer N

        let mgRecip = 1.0 / tm1000 + modA + b*l + fGC * (c+modD*l) + (1.0/(2.0 * float (N-1)) * (e+f*l+modG*l*l) )
        1.0/mgRecip |> K2C

    /// Get oligo Tm for the first N bases of oligo array
    let _temp (p:PrimerParams) (oligo: char array) (N:int) =
        // http://pubs.acs.org/doi/pdf/10.1021/bi702363u
        // figure 9
        if p.monovalentConc < 1e-8<M> then  
            // eq16 table 2
            1,999.0,eq16Table2 p.divalentConc p.primerConc p.templateConc p.dNTPConc oligo N // case 1
        else 
            let r = sqrt(p.divalentConc*1.0<M> (* UNITHACK *) ) / p.monovalentConc
            if r < 0.22 then 
                // use monovalent salt correct eq 4
                2,r,monoSaltCorrect p.monovalentConc p.primerConc p.templateConc oligo N // case 2
            else
                if r < 6.0 then 
                    // Modified eq 16 table 2
                    // with a,d,g varying according to eq 18,19,20
                    3,r,eq16Table2Mod p.divalentConc p.monovalentConc p.primerConc p.templateConc p.dNTPConc oligo N // case 3
                else
                    // eq 16 table 2
                    4,r,eq16Table2 p.divalentConc p.primerConc p.templateConc p.dNTPConc oligo N // case 4
    
    /// Get oligo Tm for the first N bases of oligo array for given parameters
    let temp (p:PrimerParams) (oligo: char array) (N:int) =
        let _,_,x = _temp p oligo N
        x


