import QtQuick 2.4
import QtQuick.Window 2.2
import QtQuick.Extras 1.4
import QtGraphicalEffects 1.0
Window { 
	title: qsTr("Grn Viewer")
	width: 200*3 + 20*3
	height: 450 
	color :"#333"
    visible: true
 	flags : Qt.WindowStaysOnTopHint | Qt.Window
 	Component {
 	    id: proteinComp 
		Rectangle {
			property real concentration 
			concentration : val < 1.25 ? val*0.8 : 1.0
			height : 20
			Text {
				anchors.left:parent.left
				anchors.leftMargin:parent.leftMargin
				text: name + ": " + val 
				color : "white"
				z : 1
			}
			width : parent.width * concentration
			anchors.left : parent.left
			color : Qt.hsla(0.4 + concentration*0.3, 0.5, 0.35, 1.0)
		}
 	}
	ListView {
		id: grnIn
		model: grnInputs
		width: 200 
		height: count * (20 + spacing)
		spacing: 5
		anchors.left: parent.left
		anchors.top: parent.top
		anchors.topMargin: 10
		anchors.leftMargin: 20
		delegate: proteinComp
	}
	ListView {
		id: grnReg
		model: grnReguls
		width: 200 
		height: count * (20 + spacing)
		spacing: 5
		anchors.left: grnIn.right
		anchors.top: parent.top
		anchors.topMargin: 10
		anchors.leftMargin: 20
		delegate: proteinComp
	}
	ListView {
		id: grnOut
		model: grnOutputs
		width: 200 
		height: count * (20 + spacing)
		spacing: 5
		anchors.left: grnReg.right
		anchors.top: parent.top
		anchors.topMargin: 10
		anchors.leftMargin: 20
		delegate: proteinComp
	}

}


