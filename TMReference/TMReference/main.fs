open System
open System.IO

let usage() =
    printfn "# Amyris reference Tm calculation implementation"
    printfn "TMReference oligoSeq  | fileOfOligosOnePerLine"
    exit 1

let procOne (s:string) =
    let s' = s.Trim().ToUpper()
    if s'.Length>0 then
        printfn "%s %A" s' (Amyris.tmcore.temp Amyris.tmcore.defaultParams (s'.ToCharArray()) s'.Length)

let validNuc (c:char) = ("ATCG"|> Set.ofSeq).Contains(c)

[<EntryPoint>]
let main argv =
    match argv with
        | [| "-h" |] | [| "--help"|] -> usage()
        | [| name |] ->
            if File.Exists(name) then
                File.ReadAllLines(name) |> Array.iter (procOne)
                0
            elif (name.ToUpper() |> Seq.forall (validNuc)) then
                procOne name
                0
            else 
                usage()
        | _ -> usage()
            