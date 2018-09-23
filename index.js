const electron = require('electron')
const {app, BrowserWindow} = electron


app.on('ready',() => {
	let win = new BrowserWindow({width:1200, height:600, icon: __dirname + '/images/ihp_icon_wc2.png'
	})
	win.loadURL(`file://${__dirname}/index.html`)
})
