
(* $Id$ *)

(*  Copyright 2004 Sascha Husa, Ian Hinder, Christiane Lechner

    This file is part of Kranc.

    Kranc is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Kranc is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kranc; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*)

BeginPackage["CodeGen`", {"Errors`", "Kranc`"}];

FlattenBlock::usage = "FlattenBlock[block] converts 'block' to a string.";
GenerateFile::usage = "GenerateFile[name, block] writes 'block' to a file of the " <>
  "specified 'name'.";
SpaceSeparated::usage = "";
NewlineSeparated::usage = "";
InfoVariable::usage = "";
CommaNewlineSeparated::usage = "";
CommaSeparated::usage = "";
Stringify::usage = "";
Quote::usage = "Quote[x] returns x surrounded by quotes";
IndentBlock::usage = "";
IndentBlock2::usage = "";
CheckBlock::usage = "";

(* In strict mode, codegen blocks must have head CodeBlock, indicating
   that they were already checked for correctness.  This helps ensure
   that blocks are checked when they are created, and errors are
   reported earlier. Kranc does not yet conform to this new style. *)

$CodeGenStrict = False;

CodeGenBlock = 
  If[$CodeGenStrict,
    CodeBlock[_],
     _CodeBlock | _String | _?AtomQ | List[(_?(MatchQ[#, CodeGenBlock] &)) ...]];

Boolean = (True | False);

CodeBlock;
CodeBlockContents;
MakeCodeBlock;

Begin["`Private`"];

(* Code generation utilities; not specific to any language *)

CheckBlock[s_String] := s;

CheckBlock[a_?AtomQ] := a;

CheckBlock[l_List] := Map[CheckBlock, l];

CheckBlock[b_CodeBlock] := CheckBlock[CodeBlockContents[b]];

FlattenBlock[b_] :=
  Module[
    {flattenBlock},
    flattenBlock[x_String] := x;
    flattenBlock[l_List] := StringJoin@@Map[FlattenBlock, l];
    flattenBlock[a_?AtomQ] := ToString[a];
    flattenBlock[x_CodeBlock] := flattenBlock[CodeBlockContents[x]];
    flattenBlock[x_] := ThrowError["Invalid arguments to flattenBlock: ", c];

    CheckBlock[b];
    flattenBlock[b]];

DefFn[
  IndentBlock[block:CodeGenBlock] :=
    StringDrop["  " <> StringReplace[FlattenBlock[block], {"\n" -> "\n  "}],-2]];

(* This should be used everywhere - need to tidy up the newline convention in CodeGen *)
DefFn[
  IndentBlock2[block:CodeGenBlock] :=
  Riffle[Map[StringJoin["  ",#] &,
      StringSplit[FlattenBlock[block],"\n"]],"\n"]];

DefFn[
  GenerateFile[filename_String, contents_] :=
  Module[
    {fp = OpenWrite[filename]},
    CheckBlock[contents];
    WriteString[fp, FlattenBlock[contents]];
    Sow[filename, GenerateFile];
    Close[fp]]];

DefFn[
  CommaNewlineSeparated[l_List] :=
  Riffle[l, ",\n"]];

DefFn[
  SpaceSeparated[l_List] :=
  Riffle[l, " "]];

DefFn[
  CommaSeparated[l_List] :=
  Riffle[l, ", "]];

DefFn[
  NewlineSeparated[l_List] :=
  Riffle[l, "\n"]];

(* Turn a section of code into a string:
   1. quote all quotes (replace all quotes with backslash-quote)
   2. break the string into lines to make it readable (replace all newlines
      with quote-newline-quote)
   3. surround the result with quotes *)
DefFn[
  Stringify[x:CodeGenBlock] :=
  "\"" <> StringReplace[StringReplace[FlattenBlock[x], "\"" -> "\\\""],
                        "\n" -> "\\n\"\n\""] <> "\"\n"];

DefFn[
  Quote[x:CodeGenBlock] :=
  {"\"", x, "\""}];

DefFn[
  MakeCodeBlock[x_] :=
  CodeBlock[CheckBlock[x]]];

DefFn[
  CodeBlockContents[CodeBlock[x_]] :=
  x];

End[];

EndPackage[];
